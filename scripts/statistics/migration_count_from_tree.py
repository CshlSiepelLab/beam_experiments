#!/usr/bin/env python3

import sys
from ete3 import Tree

def get_migration_counts(tree):
    migration_counts = {}
    for node in tree.traverse():
        if node.is_root():
            continue
        else:
            node_tissue = node.name.split("_")[-1]
            parent_tissue = node.up.name.split("_")[-1]
            if node_tissue == parent_tissue:
                continue
            migration = f"{parent_tissue}_{node_tissue}"
            if migration not in migration_counts:
                migration_counts[migration] = 1
            else:
                migration_counts[migration] += 1
    return migration_counts

# user input
true_tree_file=sys.argv[1]

# true_tree_file="results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/1983/cell_tree_seed1983_tissue_labeled_tree.nwk"

# set primary tissue
primary_tissue=sys.argv[2]

# read in tree files to ete3 tree
true_tree = Tree(true_tree_file, format=8)

# set tree root to primary
true_tree.get_tree_root().name = f'0_{primary_tissue}'

# get counts of migration events in a dict with source_recipient tissue key and count integer value
true_counts=get_migration_counts(true_tree)

# get general totals for formulas
total_true = sum(true_counts.values())

# output total migration count
print(f"Migration count: {total_true}")

