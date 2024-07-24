#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from ete3 import Tree
import dendropy
from copy import deepcopy
from arviz import hdi
import glob

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

dirs = (sys.argv[1]).split(",")
outdir = sys.argv[2]

machina_precisions = np.zeros(len(dirs))
machina_recalls = np.zeros(len(dirs))
consensus_precisions = np.zeros(len(dirs))
consensus_recalls = np.zeros(len(dirs))
random_precisions = np.zeros(len(dirs))
random_recalls = np.zeros(len(dirs))
all_thresh_rows = []
i=0

for true_tree_file in dirs:

    # get working dir from outdir based on snakemake input
    dir = os.path.dirname(os.path.dirname(os.path.dirname(true_tree_file)))
    print(f"Working dir: {dir}")

    # get sim id as dir name
    sim = os.path.basename(os.path.dirname(true_tree_file))
    print(f"Processing sim: {sim}")
    
    # get other files to compare
    machina_file = f"{dir}/machina/{sim}/machina_tree_all_tissue_labels.nwk"
    beast_posterior_file = f"{dir}/metastabayes/{sim}/combined.trees"
    consensus_file = f"{dir}/random_consensus_tissue_inference/{sim}/consensus_tissues.nwk"
    random_file = f"{dir}/random_consensus_tissue_inference/{sim}/random_tissues.nwk"
    outfile = f"{outdir}/{sim}/precision_recall.pdf"
    outfile_metrics = f"{outdir}/{sim}/metrics.csv"

    # print all file paths to output
    print(f"True tree file: {true_tree_file}")
    print(f"Machina file: {machina_file}")
    print(f"Beast posterior file: {beast_posterior_file}")
    print(f"Consensus file: {consensus_file}")
    print(f"Random file: {random_file}")

    # process true input file to get migration count dict
    if true_tree_file.endswith(".csv"):
        true_counts = process_csv(true_tree_file)
    else:
        true_counts = process_tree(true_tree_file)

    # process machina result to get precision and recall
    if os.path.exists(machina_file):
        machina_counts = process_tree(machina_file)
        machina_f1, machina_recall, machina_precision = calculate_metrics(true_counts, machina_counts)
        machina_precisions[i] = machina_precision
        machina_recalls[i] = machina_recall

    if os.path.exists(random_file):
        random_counts = process_tree(random_file)
        random_f1, random_recall, random_precision = calculate_metrics(true_counts, random_counts)
        random_precisions[i] = random_precision
        random_recalls[i] = random_recall

    if os.path.exists(consensus_file):
        consensus_counts = process_tree(consensus_file)
        consensus_f1, consensus_recall, consensus_precision = calculate_metrics(true_counts, consensus_counts)
        consensus_precisions[i] = consensus_precision
        consensus_recalls[i] = consensus_recall

    # process beast posterior result to get precision and recall
    if os.path.exists(beast_posterior_file):
        burnin_percent=0.1
        primary_tissue="P"
        beast_tree_list = dendropy.TreeList()
        beast_tree_list.read(path=beast_posterior_file, schema="nexus")
        num_beast_trees = len(beast_tree_list)
        if num_beast_trees == 0:
            continue
        num_discard = round(num_beast_trees * burnin_percent)
        beast_tree_list = beast_tree_list[num_discard:]

        posteriors = []
        # f1_scores=[]
        # precisions=[]
        # recalls=[]
        all_inferred_counts = []
        for tree in beast_tree_list:
            posterior = float(tree.annotations.get_value('posterior'))
            posteriors.append(posterior)
            remove_zero_length_nodes(tree)
            inferred_tree = dendropy_beast_to_ete_newick_with_strict_locations(tree)
            inferred_tree.get_tree_root().name = f'0_{primary_tissue}'
            inferred_counts=get_migration_counts(inferred_tree)
            # f1, recall, precision = calculate_metrics(true_counts, inferred_counts)
            # f1_scores.append(f1)
            # precisions.append(precision)
            # recalls.append(recall)
            all_inferred_counts.append(inferred_counts)

        # fit a gaussian kernel density estimate to the posterior values to get a probability density function
        pdf = gaussian_kde(posteriors)
        posterior_probs = [pdf(posterior)[0] for posterior in posteriors]
        total_posterior_prob = sum(posterior_probs)
        posterior_probs = [posterior_prob/total_posterior_prob for posterior_prob in posterior_probs]

        # calculate total counts weighted by posterior probability
        total_counts = {}
        for prob, inferred_counts in zip(posterior_probs, all_inferred_counts):
            for pattern, count in inferred_counts.items():
                for num in range(1, count+1):
                    edge = f"{pattern}_{num}"
                    if edge not in total_counts:
                        total_counts[edge] = prob
                    else:
                        total_counts[edge] += prob

        # compute thresholded precision and recall values
        thresholds = [i for i in np.arange(0, 1.00, 0.01)]
        rows = []
        for thresh in thresholds:
            thresh_counts = {key: value for key, value in total_counts.items() if value > thresh}
            edges = ['_'.join(edge.split("_")[:-1]) for edge in thresh_counts.keys()]
            thresh_counts = {}
            for edge in edges:
                if edge not in thresh_counts:
                    thresh_counts[edge] = 1
                else:
                    thresh_counts[edge] += 1
            f1, recall, precision = calculate_metrics(true_counts, thresh_counts)
            rows.append({'Threshold': thresh, 'precision': precision, 'recall': recall, 'sim': i})
        thresh_prec_rec = pd.DataFrame(rows)

    # plot precision recall curve
    size = 75
    textsize=18
    plt.figure()
    if os.path.exists(beast_posterior_file):
        plt.scatter(thresh_prec_rec['recall'], thresh_prec_rec['precision'], c=thresh_prec_rec['Threshold'], cmap='viridis', s=size, marker='x')
        plt.plot(thresh_prec_rec['recall'], thresh_prec_rec['precision'], color = "grey", label='Posterior')
    if os.path.exists(machina_file):
        plt.scatter(machina_recall, machina_precision, color='red', label='Machina', s=size, marker = "x")
    if os.path.exists(consensus_file):
        plt.scatter(consensus_recall, consensus_precision, color='blue', label='Consensus', s=size, marker = "x")
    if os.path.exists(random_file):
        plt.scatter(random_recall, random_precision, color='black', label='Random', s=size, marker = "x")
    plt.xlim(-0.05,1.05)
    plt.ylim(-0.05,1.05)
    plt.xlabel('Recall', fontsize=textsize)
    plt.ylabel('Precision', fontsize=textsize)
    plt.xticks(fontsize=textsize)
    plt.yticks(fontsize=textsize)
    plt.legend(bbox_to_anchor=(1.05, 0.4), loc='upper left', fontsize=14, edgecolor='none')
    cbar = plt.colorbar(shrink=0.4, orientation="vertical", drawedges=False, anchor=(1.05, 0.80))
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('Posterior threshold', fontsize=14, rotation=90)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


    all_thresh_rows.extend(rows)
    i+=1

    # # write metrics used for the plot to a file
    # sim
    # machina_f1, machina_recall, machina_precision
    # random_f1, random_recall, random_precision
    # consensus_f1, consensus_recall, consensus_precision


# make an overall averaged precision/recall curve
outfile = f"{outdir}/precision_recall.pdf"
avg_machina_precision = sum(machina_precisions) / len(machina_precisions)
avg_machina_recall = sum(machina_recalls) / len(machina_recalls)
avg_random_precision = sum(random_precisions) / len(random_precisions)
avg_random_recall = sum(random_recalls) / len(random_recalls)
avg_consensus_precision = sum(consensus_precisions) / len(consensus_precisions)
avg_consensus_recall = sum(consensus_recalls) / len(consensus_recalls)
all_thresh_df = pd.DataFrame(all_thresh_rows)
avg_df = all_thresh_df.groupby('Threshold')[['precision', 'recall']].mean().reset_index()

size = 75
textsize = 18
plt.figure()
plt.scatter(avg_df['recall'], avg_df['precision'], c=avg_df['Threshold'], cmap='viridis', s=size, marker='x')
plt.plot(avg_df['recall'], avg_df['precision'], color = "grey", label='Posterior')
plt.scatter(avg_machina_recall, avg_machina_precision, color='red', label='Machina', s=size, marker = "x")
plt.scatter(avg_consensus_recall, avg_consensus_precision, color='blue', label='Consensus', s=size, marker = "x")
plt.scatter(avg_random_recall, avg_random_precision, color='black', label='Random', s=size, marker = "x")
plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel('Recall', fontsize=textsize)
plt.ylabel('Precision', fontsize=textsize)
plt.xticks(fontsize=textsize)
plt.yticks(fontsize=textsize)
plt.legend(bbox_to_anchor=(1.05, 0.4), loc='upper left', fontsize=14, edgecolor='none')
cbar = plt.colorbar(shrink=0.4, orientation="vertical", drawedges=False, anchor=(1.05, 0.80))
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_ylabel('Posterior threshold', fontsize=14, rotation=90)
plt.tight_layout()
plt.savefig(outfile)
plt.close()