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

def calculate_metrics(true_counts, inferred_counts):
    TP = 0
    FP = 0
    FN = 0
    all_keys = set(true_counts.keys()).union(set(inferred_counts.keys()))
    for key in all_keys:
        if key in inferred_counts:
            inferred_count = inferred_counts[key]
            if key in true_counts:
                true_count = true_counts[key]
                if inferred_count >= true_count:
                    TP += true_count
                    FP += inferred_count - true_count
                else:
                    TP += inferred_count
                    FN += true_count - inferred_count
            else:
                FP += inferred_count
        else:
            FN += true_counts[key]
    # compute precision as TP/(TP + FP) and recall as TP/(TP + FN)
    if (TP + FP) != 0:
        precision = TP/(TP + FP)
    else:
        precision = 0
    if (TP + FN) != 0:
        recall = TP/(TP + FN)
    else:
        recall = 0
    # calculate F1 score (2((precision * recall)/(precision + recall)))
    if precision + recall == 0:
        f1 = 0
    else:
        f1 = 2 * ((precision * recall) / (precision + recall))
    return f1, recall, precision


# user input filepaths
true_tree_file=sys.argv[1] # can also be a csv of source to recipient connections to bypass the tree processing steps
inferred_tree_file=sys.argv[2]

# true_tree_file="/grid/siepel/home_norepl/staklins/stephen_data/beast_migration_inference/individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24/mS/5926/cell_tree_seed1180317166_tissue_labeled_tree.nwk"
# inferred_tree_file="/grid/siepel/home_norepl/staklins/stephen_data/beast_migration_inference/individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24/mS/5926/machina_tree_all_tissue_labels.nwk"


# process true input file to get migration count dict
if true_tree_file.endswith(".csv"):
    true_counts = process_csv(true_tree_file)
else:
    true_counts = process_tree(true_tree_file)

# process inferred input as tree to get migration count dict
inferred_counts = process_tree(inferred_tree_file)

f1, precision, recall = calculate_metrics(true_counts, inferred_counts)
print(f"F1 score: {f1} Precision: {precision} Recall: {recall}")
