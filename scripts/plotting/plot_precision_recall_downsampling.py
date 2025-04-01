#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from ete3 import Tree
import dendropy
from copy import deepcopy
from arviz import hdi
import glob
import ast
import multiprocessing as mp
from sklearn.metrics import auc


# default colors taken from metient method for consistency in visualizations
DEFAULT_COLORS = [
    "#6aa84f",
    "#be5742e1",
    "#6fa8dc",
    "#e69138",
    "#9e9e9e",
    "#c27ba0",
    "brown",
    "black",
    "darkgreen",
    "purple",
    "blue",
] * 3


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


def optionally_add_origin_migration(tree, counts, primary):
    """
    Adds the origin branch to the migration counts if the root node location is different from the known origin primary tissue.

    Parameters:
    tree (ete3.Tree): The phylogenetic tree object.
    counts (dict): A dictionary containing migration counts.
    primary (str): The known primary tissue for the origin.

    Returns:
    dict: Updated migration counts including the origin branch if applicable.
    """
    counts_copy = counts
    root_node_location = tree.get_tree_root().name.split("_")[-1]
    if root_node_location != primary:
        migration = f"{primary}_{root_node_location}"
        if migration not in counts_copy:
            counts_copy[migration] = 1
        else:
            counts_copy[migration] += 1
    return counts_copy


def process_tree(filepath):
    # set primary tissue
    primary_tissue = "P"
    # read in tree files to ete3 tree
    tree = Tree(filepath, format=8)
    # set tree root to primary
    tree.get_tree_root().name = f"0_{primary_tissue}"
    # get counts of migration events in a dict with source_recipient tissue key and count integer value
    counts = get_migration_counts(tree)
    return counts


def process_csv(filepath):
    # read in csv file to dict
    counts = {}
    with open(filepath, "r") as f:
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
            if parent is not None:  # prevents root from being lost
                parent.remove_child(node)
                children = node.child_nodes()
                for child in children:
                    parent.add_child(child)


def dendropy_beast_to_ete_newick_with_strict_locations(tree):
    # Note: this method does not convert the branch lengths properly since they are not used
    tree_copy = deepcopy(tree)
    i = 0
    for node in tree_copy.preorder_node_iter():
        try:
            prediction = node.taxon.label + "_" + node.annotations.get_value("location")
            node.taxon.label = prediction
        except Exception as e:
            prediction = f"node{i}" + "_" + node.annotations.get_value("location")
            i += 1
        node.label = prediction
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("'", ""), format=3)
    return ete_tree


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
        precision = TP / (TP + FP)
    else:
        precision = 0
    if (TP + FN) != 0:
        recall = TP / (TP + FN)
    else:
        recall = 0
    # calculate F1 score (2((precision * recall)/(precision + recall)))
    if precision + recall == 0:
        f1 = 0
    else:
        f1 = 2 * ((precision * recall) / (precision + recall))
    return f1, recall, precision


def posterior_threshold_metrics(all_inferred_counts, true_counts, i):
    # calculate total counts weighted by posterior probability
    post_prob_precision = 0
    post_prob_recall = 0
    post_prob_f1 = 0
    total_counts = {}
    prob = 1 / len(all_inferred_counts)
    for inferred_counts in all_inferred_counts:
        for pattern, count in inferred_counts.items():
            for num in range(1, count + 1):
                edge = f"{pattern}_{num}"
                if edge not in total_counts:
                    total_counts[edge] = prob
                else:
                    total_counts[edge] += prob
        # get posterior prob weighted precision, recall, and f1
        f1, recall, precision = calculate_metrics(true_counts, inferred_counts)
        post_prob_precision += prob * precision
        post_prob_recall += prob * recall
        post_prob_f1 += prob * f1

    # compute thresholded precision and recall values
    max_threshold = max(list(total_counts.values()))
    thresholds = [i for i in np.arange(0, max_threshold, 0.01)]
    rows = []
    for thresh in thresholds:
        thresh_counts = {
            key: value for key, value in total_counts.items() if value > thresh
        }
        edges = ["_".join(edge.split("_")[:-1]) for edge in thresh_counts.keys()]
        thresh_counts = {}
        for edge in edges:
            if edge not in thresh_counts:
                thresh_counts[edge] = 1
            else:
                thresh_counts[edge] += 1
        f1, recall, precision = calculate_metrics(true_counts, thresh_counts)
        rows.append(
            {
                "Threshold": thresh,
                "precision": precision,
                "recall": recall,
                "sim": i,
                "thresh_counts": thresh_counts,
            }
        )
    thresh_prec_rec = pd.DataFrame(rows)

    return thresh_prec_rec, rows, total_counts


# user inputs
dirs = (sys.argv[1]).split(",")
primary_tissue = sys.argv[2]
outdir = sys.argv[3]

# set based on simulated dataset conditions
downsampling_thresholds = [0.0, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0]
threads = 8

avg_dfs = []


def process(args):
    ds, outdir, dirs, primary_tissue = args

    all_thresh_rows = []
    i = 0

    for true_tree_file in dirs:

        # makes sure true_tree_file is not an empty string
        if not true_tree_file:
            continue

        # get working dir from outdir based on snakemake input
        dir = os.path.dirname(os.path.dirname(os.path.dirname(true_tree_file)))

        # get sim id as dir name
        sim = os.path.basename(os.path.dirname(true_tree_file))

        # get other files to compare
        beast_posterior_file = f"{dir}/metastabayes/{sim}_{ds}/combined.trees"

        outfile = f"{outdir}/{sim}/precision_recall.pdf"

        # process true input file to get migration count dict
        try:
            if true_tree_file.endswith(".csv"):
                true_counts = process_csv(true_tree_file)
            else:
                true_counts = process_tree(true_tree_file)
        except Exception as e:
            print(f"Error processing {true_tree_file}: {e}")
            continue

        # process beast posterior result to get precision and recall
        post_prob_f1 = float("nan")
        post_prob_recall = float("nan")
        post_prob_precision = float("nan")
        if os.path.exists(beast_posterior_file):
            burnin_percent = 0.1
            beast_tree_list = dendropy.TreeList()
            beast_tree_list.read(path=beast_posterior_file, schema="nexus")
            num_beast_trees = len(beast_tree_list)
            if num_beast_trees == 0:
                continue
            num_discard = round(num_beast_trees * burnin_percent)
            beast_tree_list = beast_tree_list[num_discard:]

            posteriors = []
            all_inferred_counts = []
            for tree in beast_tree_list:
                posterior = float(tree.annotations.get_value("posterior"))
                posteriors.append(posterior)
                remove_zero_length_nodes(tree)
                inferred_tree = dendropy_beast_to_ete_newick_with_strict_locations(tree)
                inferred_counts = get_migration_counts(inferred_tree)
                inferred_counts = optionally_add_origin_migration(
                    inferred_tree, inferred_counts, primary_tissue
                )
                all_inferred_counts.append(inferred_counts)

            thresh_prec_rec, rows, posterior_prob_graph = posterior_threshold_metrics(
                all_inferred_counts, true_counts, sim
            )
            all_thresh_rows.extend(rows)

        i += 1

    all_thresh_df = pd.DataFrame(all_thresh_rows)
    if not all_thresh_df.empty:
        all_thresh_df.to_csv(f"{outdir}/{ds}/all_threshold_stats.csv", index=False)
        avg_df = (
            all_thresh_df.groupby("Threshold")[["precision", "recall"]]
            .mean()
            .reset_index()
        )
    return ds, avg_df


# if we want to compute the precision recall values from scratch
# with mp.Pool(threads) as pool:
#     avg_dfs = pool.map(process, [(ds, outdir, dirs, primary_tissue) for ds in downsampling_thresholds])

# if compute is done previously, so we only want to plot precision recall from reading in the threshold csv files
filepaths = [f"{outdir}/{ds}/all_threshold_stats.csv" for ds in downsampling_thresholds]
avg_dfs = []
for filepath in filepaths:
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
        avg_df = df.groupby("Threshold")[["precision", "recall"]].mean().reset_index()
        ds = float(os.path.basename(os.path.dirname(filepath)))
        avg_dfs.append((ds, avg_df))

size = 75
textsize = 18
plt.figure()
colors = DEFAULT_COLORS
# Calculate AUC for each curve and plot those values on the plot
for i, (ds, avg_df) in enumerate(avg_dfs):
    if not avg_df.empty:
        recall = avg_df["recall"]
        precision = avg_df["precision"]
        plt.plot(recall, precision, color=colors[i], label=f"{ds}")
        # plt.scatter(avg_df['recall'], avg_df['precision'], c=avg_df['Threshold'], cmap='viridis', s=25, marker='x')
plt.xlim(-0.05, 1.05)
# plt.xlim(0.4, 1.01)
# plt.xticks(np.arange(0.4, 1.01, 0.2))
plt.ylim(-0.05, 1.05)
plt.xlabel("Recall", fontsize=textsize)
plt.ylabel("Precision", fontsize=textsize)
plt.xticks(fontsize=textsize)
plt.yticks(fontsize=textsize)
plt.legend(
    title="Downsampling\ndistance\nthreshold",
    bbox_to_anchor=(1.05, 0.8),
    loc="upper left",
    fontsize=14,
    title_fontsize=14,
    edgecolor="none",
)
# cbar = plt.colorbar(shrink=0.4, orientation="vertical", drawedges=False, anchor=(1.05, 0.80))
# cbar.ax.tick_params(labelsize=14)
# cbar.ax.set_ylabel('Posterior threshold', fontsize=14, rotation=90)
plt.tight_layout()
outfile = f"{outdir}/precision_recall.pdf"
plt.savefig(outfile)
plt.close()
