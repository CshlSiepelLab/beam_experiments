#!/usr/bin/env python3

import sys
from ete3 import Tree
import dendropy
from copy import deepcopy
import numpy as np
from arviz import hdi
import pandas as pd

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

def remove_zero_length_nodes(tree):
    for node in tree.internal_nodes():
        if node.edge_length == 0:
            parent = node.parent_node
            if parent is not None:
                parent.remove_child(node)
                children = node.child_nodes()
                for child in children:
                    parent.add_child(child)

def dendropy_beast_to_ete_newick_with_strict_locations(tree):
    tree_copy = deepcopy(tree)
    i = 0
    for node in tree_copy.preorder_node_iter():
        try:
            prediction = node.taxon.label + "_" + node.annotations.get_value('location')
            node.taxon.label = prediction
        except Exception as e:
            prediction = f"node{i}" + "_" + node.annotations.get_value('location')
            i += 1
        node.label = prediction
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("\'", ""), format=3)
    return ete_tree

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

def calculate_metrics(true_counts, inferred_counts):
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
    return f1, recall, precision


true_tree_file=sys.argv[1] # can also be csv of source,recipient format
beast_trees_file=sys.argv[2]

# true_tree_file="/grid/siepel/home_norepl/staklins/stephen_data/beast_migration_inference/individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24/mS/5926/migration_graph_seed1180317166.csv"
# beast_trees_file="/grid/siepel/home_norepl/staklins/stephen_data/beast_migration_inference/individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24/mS/5926/joint_inference_beast_combined_tissues.trees"

burnin_percent=0.1
primary_tissue="P"

# true counts are the same for all beast tree comparisons
# process true input file to get migration count dict
if true_tree_file.endswith(".csv"):
    true_counts = process_csv(true_tree_file)
else:
    true_counts = process_tree(true_tree_file)

beast_tree_list = dendropy.TreeList()
beast_tree_list.read(path=beast_trees_file, schema="nexus")

num_beast_trees = len(beast_tree_list)
num_discard = round(num_beast_trees * burnin_percent)
beast_tree_list = beast_tree_list[num_discard:]

posteriors = []
f1_scores=[]
precisions=[]
recalls=[]
rows = []
posterior_inferred_counts = {}

for tree in beast_tree_list:
    posterior = float(tree.annotations.get_value('posterior'))
    posteriors.append(posterior)
    remove_zero_length_nodes(tree)
    inferred_tree = dendropy_beast_to_ete_newick_with_strict_locations(tree)

    # calculate beast counts for each tree, get metrics, and weight by posterior probability
    inferred_tree.get_tree_root().name = f'0_{primary_tissue}'
    inferred_counts=get_migration_counts(inferred_tree)

    f1, recall, precision = calculate_metrics(true_counts, inferred_counts)
    f1_scores.append(f1)
    precisions.append(precision)
    recalls.append(recall)

    posterior_precision_recalls_row = {'Posterior': posterior, 'Precision': precision, 'Recall': recall, 'F1': f1}
    rows.append(posterior_precision_recalls_row)

    # make dict of psoterior value to dict of migration graph route counts for 95% CI calculation later
    match_true_graph = inferred_counts == true_counts
    posterior_inferred_counts[posterior] = match_true_graph

posterior_precision_recalls = pd.DataFrame(rows)

avg_posterior_f1 = sum(f1_scores) / len(f1_scores)
avg_posterior_precision = sum(precisions) / len(precisions)
avg_posterior_recall = sum(recalls) / len(recalls)
print(f"F1 score: {avg_posterior_f1} Precision: {avg_posterior_precision} Recall: {avg_posterior_recall}")

# output all precision and recall values for all posteriors
outfile = beast_trees_file.replace(".trees", "_trees_precision_recall.csv")
posterior_precision_recalls.to_csv(outfile, index=False)

### Credible interval
# # for equal tailed credible interval
# lower, upper = np.percentile(posteriors, [2.5, 97.5])
# inferred_counts_95ci = {k: v for k, v in sorted_posterior_inferred_counts.items() if lower < k < upper}

# for highest posterior density credible interval
hpd_interval = hdi(np.array(posteriors), hdi_prob=0.95)
inferred_counts_95ci = {k: v for k, v in posterior_inferred_counts.items() if hpd_interval[0] <= k <= hpd_interval[1]}

num_perfect_graphs = sum(inferred_counts_95ci.values())
print(f"Number of perfect graphs in 95% HPDI CI: {num_perfect_graphs}")

found_perfect_in_95ci = False
if num_perfect_graphs > 0:
    found_perfect_in_95ci = True
print(f"Found a perfect graph in 95% HPDI CI: {found_perfect_in_95ci}")

