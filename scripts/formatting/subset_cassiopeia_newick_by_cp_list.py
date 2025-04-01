#!/usr/bin/env python3

### This script is meant to take in a newick string of many taxa and then prune the tree to result in a new newick string with only the taxa specified by another input .txt file with taxa names on newlines to keep from the original newick string

from ete3 import Tree


def prune_tree(tree, taxa_list):
    for leaf in tree.iter_leaves():
        if leaf.name not in taxa_list:
            leaf.detach()


# Replace 'newick_string' and 'taxa_list' with your actual Newick string and taxa list
cassiopeia_file = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/raw_data/MMUS1495_cassiopeia_greedy.newick"
taxa_list_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/results/cp01/cp01_asv_names.txt"

# Parse the Newick string into an ETE Tree object
tree = Tree(cassiopeia_file)

with open(taxa_list_path, "r") as file:
    taxa_list = [line.strip() for line in file.readlines()]

# Prune the tree to keep only the specified taxa
# prune_tree(tree, taxa_list)
mrca = tree.get_common_ancestor(taxa_list)


# Get the Newick string for the pruned tree
output_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/results/cp01/cp01_MMUS1495_cassiopeia_greedy.newick"
# tree.write(outfile=output_path, format=9)
mrca.write(outfile=output_path, format=9)
