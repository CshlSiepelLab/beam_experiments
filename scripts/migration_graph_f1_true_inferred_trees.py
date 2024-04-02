#!/usr/bin/env python3

import sys
from ete3 import Tree

def get_migration_counts(tree):
    migration_counts = {}
    for node in tree.traverse():
        if node.is_root() or node.is_leaf():
            continue
        else:
            node_tissue = node.name.split("_")[1]
            parent_tissue = node.up.name.split("_")[1]
            if node_tissue == parent_tissue:
                continue
            migration = f"{parent_tissue}_{node_tissue}"
            if migration not in migration_counts:
                migration_counts[migration] = 1
            else:
                migration_counts[migration] += 1
    return migration_counts

# true_tree_file=sys.argv[1]
# inferred_tree_file=sys.argv[2]

true_tree_file=""
inferred_tree_file=""

# set primary tissue
primary_tissue="P"

# read in tree files to ete3 tree
true_tree = Tree(true_tree_file, format=8)
inferred_tree = Tree(inferred_tree_file, format=8)

# set tree root to primary
true_tree.get_tree_root().name = f'0_{primary_tissue}'
inferred_tree.get_tree_root().name = f'0_{primary_tissue}'

# get counts of migration events in a dict with source_recipient tissue key and count integer value
true_counts=get_migration_counts(true_tree)
inferred_counts=get_migration_counts(inferred_tree)

# calculate accuracy

# calculate precision

# calculate F1 score
