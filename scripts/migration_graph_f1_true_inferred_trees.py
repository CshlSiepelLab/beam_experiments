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

true_tree_file=sys.argv[1]
inferred_tree_file=sys.argv[2]

# true_tree_file="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874_tissue_labeled_tree.nwk"
# inferred_tree_file="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/machina_tree_all_tissue_labels.nwk"

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

# get general totals for formulas
total_inferred = sum(inferred_counts.values())
total_true = sum(true_counts.values())

# get union of positives for true and inferred
union_positives = 0
for key in inferred_counts.keys():
    inferred_value = inferred_counts[key]
    if key not in true_counts:
        continue
    else:
        true_value = true_counts[key]
        if true_value < inferred_value:
            union_positives += true_value
        else:
            union_positives += inferred_value

# calculate recall (union of inferred and true positives over inferred total positives)
recall = union_positives / total_inferred

# calculate precision (union of inferred and true positives over total true positives)
precision = union_positives / total_true

# calculate F1 score (2((precision * recall)/(precision + recall)))
f1 = 2 * ((precision * recall) / (precision + recall))

print(f"F1 score: {f1}")
