#!/usr/bin/env python3

from ete3 import Tree
import sys

newick_file = sys.argv[1]
tips_to_keep_csv = sys.argv[2]
output_file = sys.argv[3]

# Read the tree from the newick file
tree = Tree(newick_file, format=5)

# Parse the CSV string of tips to keep
tips_to_keep = tips_to_keep_csv.split(',')

# Iterate over all leaves and prune those not in the tips_to_keep list

tree = tree.prune(tips_to_keep)

# Name nodes
i = 0
for node in tree.traverse():
    if not node.is_leaf():
        node.name = f"node{1}"
        i += 1

# Write the pruned tree to the output file
tree.write(outfile=output_file, format=8)