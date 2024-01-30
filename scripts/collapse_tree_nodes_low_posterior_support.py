#!usr/bin/env python3

# This script takes in a newick stirng tree with posterior support values for all nodes and then collapses nodes with low support to a single node.

import sys
import ete3

nexus_file = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/results/cp01/cp01_mutation_matrix_reformatted_tidetree_tidetree_sequences_formatted_for_tidetree.1706295631518.tree"

tree = ete3.Tree(nexus_file)
tree = ete3.Nexml()
tree.build_from_file(nexus_file)
