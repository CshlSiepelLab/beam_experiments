import sys
import dendropy
import numpy as np


def read_tree(tree_file):
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=tree_file, schema="newick", preserve_underscores=True, taxon_namespace=tns)
    return tree


def reset_taxon_namespace(tree, tns):
    newick_str = tree.as_string(schema="newick")
    tree = dendropy.Tree.get(data=newick_str, schema="newick", taxon_namespace=tns)
    return tree


def preprocess_tree(tree, tns):
    for node in tree:
        if not node.is_leaf():
            node.label = None
        # Just doing this for backwards compatibility with older simulated trees that have appended "_sampled" to leaf names, which we want to ignore here
        else:
            if node.label is not None and "_sampled" in node.label:
                node.label = node.label.replace("_sampled", "")
            if node.taxon.label is not None and "_sampled" in node.taxon.label:
                node.taxon.label = node.taxon.label.replace("_sampled", "")
    tree = reset_taxon_namespace(tree, tns)     
    tree.resolve_polytomies()
    tree.suppress_unifurcations()
    tree.is_rooted = True
    tree.encode_bipartitions()
    return tree


def preprocess_tree_file(tree_file, tns=None):
    tree = read_tree(tree_file)
    if tns is None:
        tns = tree.taxon_namespace
    tree = preprocess_tree(tree, tns)
    return tree


def calculate_rf_distance(tree1, tree2):
    ref_bipartitions = set(tree1.bipartition_encoding)
    comparison_bipartitions = set(tree2.bipartition_encoding)
    is_true_tree = ref_bipartitions == comparison_bipartitions
    # Calculate the Robinson-Foulds distance as number of bipartitions that differ between trees, so FP + TP from comparions tree1 bipartitions as refereence and tree3 bipartitions as comparison)
    # This calculation is from the dendropy treecompare symmetric distance function, implemented here natively
    rf_distance = len(ref_bipartitions.symmetric_difference(comparison_bipartitions))
    # Calculate normalized RF distance (between 0 and 1) with the total number of bipartitions in both trees as denominator
    nontrivial_ref_bipartitions = [b for b in ref_bipartitions if not b.is_trivial()]
    nontrivial_comparison_bipartitions = [b for b in comparison_bipartitions if not b.is_trivial()]
    num_nontrivial_bipartitions = len(nontrivial_ref_bipartitions) + len(nontrivial_comparison_bipartitions)
    normalized_rf = rf_distance / num_nontrivial_bipartitions
    return rf_distance, normalized_rf, is_true_tree


true_tree_file = sys.argv[1]    # Simulated nwk file from direct output
laml_tree_file = sys.argv[2]     # LAML inferred tree in nwk format from direct output
cassiopeia_tree_file = sys.argv[3]  # Cassiopeia-greedy inferred tree in nwk format
beam_trees_file = sys.argv[4]   # BEAM posterior trees in nexus format
outfile = sys.argv[5]
pathfinder_tree_file = sys.argv[6]  # PathFinder nj inferred tree in nwk format (optional)

true_tree = preprocess_tree_file(true_tree_file)
laml_tree = preprocess_tree_file(laml_tree_file, tns=true_tree.taxon_namespace)
cassiopeia_tree = preprocess_tree_file(cassiopeia_tree_file, tns=true_tree.taxon_namespace)

# Pathfinder tree is optional since it's not always run and also has an added Normal node that needs to be handled
if pathfinder_tree_file.lower() != "none":
    path_tns = dendropy.TaxonNamespace()
    pathfinder_tree = dendropy.Tree.get(path=pathfinder_tree_file, schema="newick", preserve_underscores=True, taxon_namespace=path_tns)
    normal_node = pathfinder_tree.find_node_with_taxon_label("Normal")
    if normal_node is None:
        raise ValueError("Node labeled 'Normal' not found in tree.")
    pathfinder_tree.reroot_at_node(normal_node)
    pathfinder_tree.prune_taxa_with_labels(["Normal"])
    pathfinder_tree = preprocess_tree(pathfinder_tree, true_tree.taxon_namespace)
    pathfinder_rf, pathfinder_normalized_rf, pathfinder_is_true_tree = calculate_rf_distance(true_tree, pathfinder_tree)
else:
    pathfinder_rf, pathfinder_normalized_rf, pathfinder_is_true_tree = None, None, None

# true_tree_rf, true_tree_normalized_rf, true_tree_is_true_tree = calculate_rf_distance(true_tree, true_tree) # Just using as a control for the RF calculation function
laml_rf, laml_normalized_rf, laml_is_true_tree = calculate_rf_distance(true_tree, laml_tree)
cassiopeia_rf, cassiopeia_normalized_rf, cassiopeia_is_true_tree = calculate_rf_distance(true_tree, cassiopeia_tree)

tns = dendropy.TaxonNamespace()
beam_trees = dendropy.TreeList.get(path=beam_trees_file, schema="nexus", preserve_underscores=True, taxon_namespace=tns)

burnin_proportion = 0.10
num_input_trees = len(beam_trees)
burnin = int(burnin_proportion * num_input_trees)
beam_trees = beam_trees[burnin:]
num_remaining_trees = len(beam_trees)
per_beam_tree_prob = 1/num_remaining_trees

posterior_rf = 0
posterior_normalized_rf = 0
beam_contains_true_tree = False
for beam_tree in beam_trees:
    beam_tree = preprocess_tree(beam_tree, true_tree.taxon_namespace)
    beam_rf, beam_normalized_rf, beam_is_true_tree = calculate_rf_distance(true_tree, beam_tree)
    posterior_rf += beam_rf * per_beam_tree_prob
    posterior_normalized_rf += beam_normalized_rf * per_beam_tree_prob
    if beam_is_true_tree:
        beam_contains_true_tree = True


# Create a random binary rooted control tree with same taxa as a control measure
num_tips = len(true_tree.leaf_nodes())
random_tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=num_tips, taxon_namespace=true_tree.taxon_namespace)
random_tree.randomly_assign_taxa()
random_tree = preprocess_tree(random_tree, true_tree.taxon_namespace)
random_rf, random_normalized_rf, random_is_true_tree = calculate_rf_distance(true_tree, random_tree)

with open(outfile, 'w') as f:
    f.write("laml_rf,laml_normalized_rf,laml_is_true_tree,cassiopeia_rf,cassiopeia_normalized_rf,cassiopeia_is_true_tree,beam_posterior_rf,beam_posterior_normalized_rf,beam_contains_true_tree,pathfinder_rf,pathfinder_normalized_rf,pathfinder_is_true_tree,random_rf,random_normalized_rf,random_is_true_tree\n")
    f.write(f"{laml_rf},{laml_normalized_rf},{laml_is_true_tree},{cassiopeia_rf},{cassiopeia_normalized_rf},{cassiopeia_is_true_tree},{posterior_rf},{posterior_normalized_rf},{beam_contains_true_tree},{pathfinder_rf},{pathfinder_normalized_rf},{pathfinder_is_true_tree},{random_rf},{random_normalized_rf},{random_is_true_tree}\n")

    