#!usr/bin/env python3

# This script takes in a newick stirng tree with posterior support values for all nodes and then collapses nodes with low support to a single node.

import sys
import dendropy

def collapse_low_support_nodes(tree, threshold):
    labels = []
    i = 0
    j = 0
    for node in tree.internal_nodes():
        j = j + 1
        posterior = node.annotations.get_value('posterior')
        if posterior is not None and float(posterior) < threshold:
            node.collapse_clade()
            i = i + 1
    print(f"Collapsed {i} out of {j} internal nodes with low posterior support.")

nexus_tree_file_path = sys.argv[1]
threshold = float(sys.argv[2])

#nexus_tree_file_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/results/cp01/cp01_mutation_matrix_reformatted_tidetree_tidetree_sequences_formatted_for_tidetree.1706295631518.tree"


tree = dendropy.Tree.get(path=nexus_tree_file_path, schema='nexus')
print("Input tree topology:")
print(tree.as_ascii_plot())
#collapse_low_support_nodes(tree, 0.5)
collapse_low_support_nodes(tree, threshold)

# Print the modified tree
print("Collapsed tree topology:")
print(tree.as_ascii_plot())
#print(tree.as_string('newick'))
output_nexus_file_path = ".".join(nexus_tree_file_path.split(".")[:-1]) + f"_collapsed_nodes.tree"
tree.write_to_path(output_nexus_file_path, schema='nexus')
