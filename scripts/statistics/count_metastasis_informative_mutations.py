import sys
import pandas as pd
from ete3 import Tree


# tree_file = sys.argv[1] # nwk file
# tissues_file = sys.argv[2]  # vertex.labeling file
# indel_matrix_file = sys.argv[3]  # indel matrix tsv file
# output_file = sys.argv[4]  # output file to write results to

tree_file = "/grid/siepel/home/staklins/projects/crispr_barcode/data/variable_migration_and_mutation_rates_8_19_24/mig4_mut0025_2528/cell_tree_seed2528.nwk"
tissues_file = "/grid/siepel/home/staklins/projects/crispr_barcode/data/variable_migration_and_mutation_rates_8_19_24/mig4_mut0025_2528/cell_tree_seed2528.vertex.labeling"
indel_matrix_file = "/grid/siepel/home/staklins/projects/crispr_barcode/data/variable_migration_and_mutation_rates_8_19_24/mig4_mut0025_2528/mig4_mut0025_2528_indel_character_matrix.tsv"
output_file = "/grid/siepel/home/staklins/projects/crispr_barcode/data/variable_migration_and_mutation_rates_8_19_24/mig4_mut0025_2528/cell_tree_seed2528_mig_informative_mutation_counts.tsv"

tree = Tree(tree_file, format=1)
for node in tree.traverse():
    node.name = str(node.name)

tissues = {}
with open(tissues_file) as f:
    for line in f:
        node, tissue = line.strip().split()
        tissues[node] = tissue

indels = pd.read_csv(indel_matrix_file, sep="\t", index_col=0).astype(int)
indels.index = indels.index.astype(str)

# Infer internal node mutation states by leveraging the irreversibility of indels and ignoring missing data as an approximation of where mutations occurred
# First do a postorder traversal to get the intersection of all child mutation sets, effectively getting the accumulation of mutations from the root to the tips
states = {}
for node in tree.traverse("postorder"):
    if node.up is None:
        states[str(node.name)] = set()  # My simulated data has an origin node with no branch length above it, so this root read in to the ete3 tree cannot have mutations
        continue
    name = str(node.name)
    if node.is_leaf():
        states[name] = {f"{col}:{mut}" for col, mut in indels.loc[name].items() if mut not in (0, -1)}
    else:
        children = node.children
        states[name] = set(states[str(children[0].name)])
        for child in children[1:]:
            states[name] = states[name].intersection(states[str(child.name)])

# Now do a preorder traversal to prune any mutations that could not have occurred on the branch leading to this node
pruned_states = {}
for node in tree.traverse("preorder"):
    name = str(node.name)
    parent = node.up
    if parent is not None:
        parent_name = str(parent.name)
        # Only keep mutations not present in the parent
        pruned_states[name] = states[name] - states[parent_name]
    else:
        pruned_states[name] = states[name]

# Now find which nodes are recipients of migration events
migration_recipients = []
for node in tree.traverse():
    parent = node.up
    if parent is not None:
        parent_name = str(parent.name)
        child_name = str(node.name)
        if tissues[parent_name] != tissues[child_name]:
            migration_recipients.append(child_name)

# Now count mutations on branches leading to nodes that are migration recipients
num_mig_informative_muts = sum([len(pruned_states[child]) for child in migration_recipients])
num_total_muts = sum([len(list(muts)) for muts in pruned_states.values()])
num_migration_edges = len(migration_recipients)
total_num_edges = len(list(tree.traverse())) - 1  # exclude root edge
num_tips = len(tree.get_leaves())

with open(output_file, "w") as f:
    f.write("num_mig_informative_muts\tnum_total_muts\tnum_migration_edges\ttotal_num_edges\tnum_tips\n")
    f.write(f"{num_mig_informative_muts}\t{num_total_muts}\t{num_migration_edges}\t{total_num_edges}\t{num_tips}\n")

