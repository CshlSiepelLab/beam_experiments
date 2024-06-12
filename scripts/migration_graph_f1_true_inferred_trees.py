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

def process_tree(filepath):
    # set primary tissue
    primary_tissue="P"
    # read in tree files to ete3 tree
    tree = Tree(filepath, format=8)
    # set tree root to primary
    tree.get_tree_root().name = f'0_{primary_tissue}'
    # get counts of migration events in a dict with source_recipient tissue key and count integer value
    counts=get_migration_counts(tree)
    return counts

def process_csv(filepath):
    # read in csv file to dict
    counts = {}
    with open(filepath, 'r') as f:
        for line in f:
            # skip header line that typically is "source,recipient"
            if "source" in line:
                continue
            source, recipient = line.strip().split(",")
            migration = f"{source}_{recipient}"
            if migration not in counts:
                counts[migration] = 1
            else:
                counts[migration] += 1
    return counts

# user input filepaths
true_tree_file=sys.argv[1] # can also be a csv of source to recipient connections to bypass the tree processing steps
inferred_tree_file=sys.argv[2]

# true_tree_file="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/new_simulator_uniformTransitionProbs_6_6_24/mS/6016/migration_graph_seed1913716328.csv"
# inferred_tree_file="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/machina_tree_all_tissue_labels.nwk"


# process true input file to get migration count dict
if true_tree_file.endswith(".csv"):
    true_counts = process_csv(true_tree_file)
else:
    true_counts = process_tree(true_tree_file)

# process inferred input as tree to get migration count dict
inferred_counts = process_tree(inferred_tree_file)

# compute statistics from the migration counts dicts
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
if precision + recall == 0:
    f1 = 0
else:
    f1 = 2 * ((precision * recall) / (precision + recall))

print(f"F1 score: {f1}")
