#!/usr/bin/env python3

import sys
import os
import re
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
import ast
import pdb

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

def posterior_threshold_metrics(posterior_probs, all_inferred_counts, true_counts, i):
    # calculate total counts weighted by posterior probability
    post_prob_precision = 0
    post_prob_recall = 0
    post_prob_f1 = 0
    total_counts = {}
    for prob, inferred_counts in zip(posterior_probs, all_inferred_counts):
        for pattern, count in inferred_counts.items():
            for num in range(1, count+1):
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
    max_threshold = max(map(float, total_counts.values()))
    thresholds = [i for i in np.arange(0, max_threshold, 0.01)]
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
        rows.append({'Threshold': thresh, 'precision': precision, 'recall': recall, 'sim': i, 'thresh_counts': thresh_counts})
    thresh_prec_rec = pd.DataFrame(rows)

    # use the 70% posterior thresholded values as the final precision, recall, and f1 instead of the full posterior weighted values
    t = 0.70
    # check that the max for the data is not below the threshold for the consensus graph
    if max_threshold < t:
        t = np.floor(max_threshold * 100) / 100
    threshold_df = thresh_prec_rec[np.isclose(thresh_prec_rec['Threshold'], t)]
    post_prob_precision = threshold_df['precision'].values[0]
    post_prob_recall = threshold_df['recall'].values[0]
    post_prob_f1 = 2 * ((post_prob_precision * post_prob_recall) / (post_prob_precision + post_prob_recall))

    return thresh_prec_rec, rows, total_counts, post_prob_f1, post_prob_recall, post_prob_precision

# user inputs
dirs = (sys.argv[1]).split(",")
primary_tissue=sys.argv[2]
outdir = sys.argv[3]
plot_independent_plots = False
debug = False

# make a file to record performance statistics for all sim datasets
outfile_metrics = f"{outdir}/metrics.csv"
with open(outfile_metrics, "w") as file:
    header="sim,random_f1,random_recall,random_precision,consensus_f1,consensus_recall,consensus_precision,parsimony_f1,parsimony_recall,parsimony_precision,machina_f1,machina_recall,machina_precision,metient_f1,metient_recall,metient_precision,pathfinder_f1,pathfinder_recall,pathfinder_precision,beast_f1,beast_recall,beast_precision\n"
    file.write(header)

machina_precisions = np.zeros(len(dirs))
machina_recalls = np.zeros(len(dirs))
metient_precisions = np.zeros(len(dirs))
metient_recalls = np.zeros(len(dirs))
consensus_precisions = np.zeros(len(dirs))
consensus_recalls = np.zeros(len(dirs))
random_precisions = np.zeros(len(dirs))
random_recalls = np.zeros(len(dirs))
parsimony_precisions = np.zeros(len(dirs))
parsimony_recalls = np.zeros(len(dirs))
all_thresh_rows = []
pathfinder_all_thresh_rows = []
i=0

for true_tree_file in dirs:

    # makes sure true_tree_file is not an empty string
    if not true_tree_file:
        continue

    # get working dir from outdir based on snakemake input
    dir = os.path.dirname(os.path.dirname(os.path.dirname(true_tree_file)))
    if debug:
        print(f"Working dir: {dir}")

    # get sim id as dir name
    sim = os.path.basename(os.path.dirname(true_tree_file))
    if debug:
        print(f"Processing sim: {sim}")
    
    # get other files to compare
    machina_file = f"{dir}/machina/{sim}/machina_tree_all_tissue_labels.nwk"
    metient_file = f"{dir}/metient/{sim}/{sim}_{primary_tissue}_migration_graphs.txt"
    beast_posterior_file = f"{dir}/metastabayes/{sim}/combined.trees"
    pathfinder_posterior_file = f"{dir}/pathfinder/{sim}/clone_aln_all_output_counts.txt"
    consensus_file = f"{dir}/random_consensus_parsimony_tissue_inference/{sim}/consensus_tissues.nwk"
    random_file = f"{dir}/random_consensus_parsimony_tissue_inference/{sim}/random_tissues.nwk"
    parsimony_file = f"{dir}/random_consensus_parsimony_tissue_inference/{sim}/parsimony_tissues.nwk"
    outfile = f"{outdir}/{sim}/precision_recall.pdf"

    # print all file paths to output
    if debug:
        print(f"True tree file: {true_tree_file}")
        print(f"Machina file: {machina_file}")
        print(f"Metient file: {metient_file}")
        print(f"Beast posterior file: {beast_posterior_file}")
        print(f"Consensus file: {consensus_file}")
        print(f"Random file: {random_file}")

    # process true input file to get migration count dict
    try:
        if true_tree_file.endswith(".csv"):
            true_counts = process_csv(true_tree_file)
        else:
            true_counts = process_tree(true_tree_file)
    except Exception as e:
        print(f"Error processing {true_tree_file}: {e}")
        continue

    # output true graph to a file
    with open(f"{outdir}/{sim}_true_graph.csv", "w") as file:
        file.write(f"source_target,num_edges\n")
        for key, value in true_counts.items():
            file.write(f"{key},{value}\n")

    # process machina result to get precision and recall
    machina_f1 = float('nan')
    machina_recall = float('nan')
    machina_precision = float('nan')
    if os.path.exists(machina_file):
        machina_counts = process_tree(machina_file)
        machina_f1, machina_recall, machina_precision = calculate_metrics(true_counts, machina_counts)
        machina_precisions[i] = machina_precision
        machina_recalls[i] = machina_recall

    # process metient result to get precision and recall
    metient_f1 = float('nan')
    metient_recall = float('nan')
    metient_precision = float('nan')
    if os.path.exists(metient_file):
        with open(metient_file, 'r') as file:
            lines = file.readlines()
        top_loss_metient = lines[1].strip().split('\t')
        metient_counts_input = ast.literal_eval(top_loss_metient[1])
        metient_counts = {}
        for outer_key, inner_dict in metient_counts_input.items():
            for inner_key, value in inner_dict.items():
                if value != 0.0:
                    metient_counts[f"{outer_key}_{inner_key}"] = value
        metient_f1, metient_recall, metient_precision = calculate_metrics(true_counts, metient_counts)
        metient_precisions[i] = metient_precision
        metient_recalls[i] = metient_recall

    # process pathfinder result to get precision and recall
    pathfinder_post_prob_f1 = float('nan')
    pathfinder_post_prob_recall = float('nan')
    pathfinder_post_prob_precision = float('nan')
    if os.path.exists(pathfinder_posterior_file):
        pathfinder_raw_output = pd.read_csv(pathfinder_posterior_file, sep='\t')
        pathfinder_raw_output = pathfinder_raw_output.drop_duplicates()
        raw_posterior_probs = pathfinder_raw_output['probability'].tolist()
        # pathfinder probabilities do not sum to 1, so normalize them across all output migration graphs
        pathfinder_posterior_probs = [posterior_prob / sum(raw_posterior_probs) for posterior_prob in raw_posterior_probs]
        pathfinder_all_inferred_counts = []
        for raw_graph in pathfinder_raw_output['paths'].tolist():
            graph = raw_graph.split(';')
            cleaned_graph = {}
            for item in graph:
                item = re.sub(r'\[.*?\]', '', item)
                item = item.replace('->', '_')
                if item not in cleaned_graph:
                    cleaned_graph[item] = 1
                else:
                    cleaned_graph[item] += 1
            pathfinder_all_inferred_counts.append(cleaned_graph)

        pathfinder_thresh_prec_rec, pathfinder_rows, pathfinder_posterior_prob_graph, pathfinder_post_prob_f1, pathfinder_post_prob_recall, pathfinder_post_prob_precision = posterior_threshold_metrics(pathfinder_posterior_probs, pathfinder_all_inferred_counts, true_counts, sim)
        pathfinder_all_thresh_rows.extend(pathfinder_rows)

    # process random result to get precision and recall
    random_f1 = float('nan')
    random_recall = float('nan')
    random_precision = float('nan')
    if os.path.exists(random_file):
        random_counts = process_tree(random_file)
        random_f1, random_recall, random_precision = calculate_metrics(true_counts, random_counts)
        random_precisions[i] = random_precision
        random_recalls[i] = random_recall


    # process consensus result to get precision and recall
    consensus_f1 = float('nan')
    consensus_recall = float('nan')
    consensus_precision = float('nan')
    if os.path.exists(consensus_file):
        consensus_counts = process_tree(consensus_file)
        consensus_f1, consensus_recall, consensus_precision = calculate_metrics(true_counts, consensus_counts)
        consensus_precisions[i] = consensus_precision
        consensus_recalls[i] = consensus_recall
    
    # process greedy fitch parsimony result to get precision and recall
    parsimony_f1 = float('nan')
    parsimony_recall = float('nan')
    parsimony_precision = float('nan')
    if os.path.exists(parsimony_file):
        parsimony_counts = process_tree(parsimony_file)
        parsimony_f1, parsimony_recall, parsimony_precision = calculate_metrics(true_counts, parsimony_counts)
        parsimony_precisions[i] = parsimony_precision
        parsimony_recalls[i] = parsimony_recall

    # process beast posterior result to get precision and recall
    post_prob_f1 = float('nan')
    post_prob_recall = float('nan')
    post_prob_precision = float('nan')
    if os.path.exists(beast_posterior_file):
        burnin_percent=0.1
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
            inferred_counts = get_migration_counts(inferred_tree)
            inferred_counts = optionally_add_origin_migration(inferred_tree, inferred_counts, primary_tissue)
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
        thresh_prec_rec, rows, posterior_prob_graph, post_prob_f1, post_prob_recall, post_prob_precision = posterior_threshold_metrics(posterior_probs, all_inferred_counts, true_counts, sim)
        # output posterior_prob_graph to a file
        with open(f"{outdir}/{sim}_posterior_prob_graph.csv", "w") as file:
            posterior_prob_graph = dict(sorted(posterior_prob_graph.items(), key=lambda x: x[1], reverse=True))
            for key, value in posterior_prob_graph.items():
                file.write(f"{key},{value}\n")
        all_thresh_rows.extend(rows)

    if plot_independent_plots:
        # plot precision recall curve
        size = 75
        textsize=18
        plt.figure()
        if os.path.exists(beast_posterior_file):
            # plt.scatter(thresh_prec_rec['recall'], thresh_prec_rec['precision'], c=thresh_prec_rec['Threshold'], cmap='viridis', s=25, marker='x')
            plt.plot(thresh_prec_rec['recall'], thresh_prec_rec['precision'], color = "grey", label='Beast')
        if os.path.exists(beast_posterior_file):
            # plt.scatter(pathfinder_thresh_prec_rec['recall'], pathfinder_thresh_prec_rec['precision'], c=pathfinder_thresh_prec_rec['Threshold'], cmap='viridis', s=25, marker='x')
            plt.plot(pathfinder_thresh_prec_rec['recall'], pathfinder_thresh_prec_rec['precision'], color = "brown", label='PathFinder')
        if os.path.exists(machina_file):
            plt.scatter(machina_recall, machina_precision, color='red', label='Machina', s=size, marker = "x")
        if os.path.exists(metient_file):
            plt.scatter(metient_recall, metient_precision, color='green', label='Metient', s=size, marker = "x")
        if os.path.exists(consensus_file):
            plt.scatter(consensus_recall, consensus_precision, color='blue', label='Consensus', s=size, marker = "x")
        if os.path.exists(random_file):
            plt.scatter(random_recall, random_precision, color='black', label='Random', s=size, marker = "x")
        if os.path.exists(parsimony_file):
            plt.scatter(parsimony_recall, parsimony_precision, color='black', label='Random', s=size, marker = "x")
        plt.xlim(-0.05,1.05)
        plt.ylim(-0.05,1.05)
        plt.xlabel('Recall', fontsize=textsize)
        plt.ylabel('Precision', fontsize=textsize)
        plt.xticks(fontsize=textsize)
        plt.yticks(fontsize=textsize)
        plt.legend(bbox_to_anchor=(1.05, 0.4), loc='upper left', fontsize=14, edgecolor='none')
        # cbar = plt.colorbar(shrink=0.4, orientation="vertical", drawedges=False, anchor=(1.05, 0.80))
        # cbar.ax.tick_params(labelsize=14)
        # cbar.ax.set_ylabel('Posterior threshold', fontsize=14, rotation=90)
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

    i+=1

    # write metrics used for the plot to a file
    with open(outfile_metrics, "a") as file:
        data = f"{sim},{random_f1},{random_recall},{random_precision},{consensus_f1},{consensus_recall},{parsimony_f1},{parsimony_recall},{parsimony_precision},{machina_f1},{machina_recall},{machina_precision},{metient_f1},{metient_recall},{metient_precision},{pathfinder_post_prob_f1},{pathfinder_post_prob_recall},{pathfinder_post_prob_precision},{post_prob_f1},{post_prob_recall},{post_prob_precision}\n"
        file.write(data)

# make an overall averaged precision/recall curve
outfile = f"{outdir}/precision_recall.pdf"
avg_machina_precision = float('nan')
avg_machina_recall = float('nan')
avg_metient_precision = float('nan')
avg_metient_recall = float('nan')
avg_random_precision = float('nan')
avg_random_recall = float('nan')
avg_consensus_precision = float('nan')
avg_consensus_recall = float('nan')


if np.any(machina_precisions):
    avg_machina_precision = sum(machina_precisions) / len(machina_precisions)
if np.any(machina_recalls):
    avg_machina_recall = sum(machina_recalls) / len(machina_recalls)
if np.any(metient_precisions):
    avg_metient_precision = sum(metient_precisions) / len(metient_precisions)
if np.any(metient_recalls):
    avg_metient_recall = sum(metient_recalls) / len(metient_recalls)
if np.any(random_precisions):
    avg_random_precision = sum(random_precisions) / len(random_precisions)
if np.any(random_recalls):
    avg_random_recall = sum(random_recalls) / len(random_recalls)
if np.any(consensus_precisions):
    avg_consensus_precision = sum(consensus_precisions) / len(consensus_precisions)
if np.any(consensus_recalls):
    avg_consensus_recall = sum(consensus_recalls) / len(consensus_recalls)
if np.any(parsimony_precisions):
    avg_parsimony_precision = sum(parsimony_precisions) / len(parsimony_precisions)
if np.any(parsimony_recalls):
    avg_parsimony_recall = sum(parsimony_recalls) / len(parsimony_recalls)
all_thresh_df = pd.DataFrame(all_thresh_rows)
if not all_thresh_df.empty:
    avg_df = all_thresh_df.groupby('Threshold')[['precision', 'recall']].mean().reset_index()
    all_thresh_df.to_csv(f"{outdir}/all_threshold_stats.csv", index=False)

pathfinder_all_thresh_df = pd.DataFrame(pathfinder_all_thresh_rows)
if not pathfinder_all_thresh_df.empty:
    pathfinder_avg_df = pathfinder_all_thresh_df.groupby('Threshold')[['precision', 'recall']].mean().reset_index()

size = 75
textsize = 18
plt.figure()
if not avg_df.empty:
    # plt.scatter(avg_df['recall'], avg_df['precision'], c=avg_df['Threshold'], cmap='viridis', s=25, marker='x')
    plt.plot(avg_df['recall'], avg_df['precision'], color = 'grey', label='Beast')
if not pathfinder_avg_df.empty:
    # plt.scatter(pathfinder_avg_df['recall'], pathfinder_avg_df['precision'], c=pathfinder_avg_df['Threshold'], cmap='viridis', s=25, marker='x')
    plt.plot(pathfinder_avg_df['recall'], pathfinder_avg_df['precision'], color='brown', label='Pathfinder')
if not np.isnan(avg_machina_recall) and not np.isnan(avg_machina_precision):
    plt.scatter(avg_machina_recall, avg_machina_precision, color='red', label='Machina', s=size, marker="x")
if not np.isnan(avg_metient_recall) and not np.isnan(avg_metient_precision):
    plt.scatter(avg_metient_recall, avg_metient_precision, color='green', label='Metient', s=size, marker="x")
if not np.isnan(avg_consensus_recall) and not np.isnan(avg_consensus_precision):
    plt.scatter(avg_consensus_recall, avg_consensus_precision, color='blue', label='Consensus', s=size, marker="x")
if not np.isnan(avg_random_recall) and not np.isnan(avg_random_precision):
    plt.scatter(avg_random_recall, avg_random_precision, color='black', label='Random', s=size, marker="x")
if not np.isnan(avg_parsimony_recall) and not np.isnan(avg_parsimony_precision):
    plt.scatter(avg_parsimony_recall, avg_parsimony_precision, color='purple', label='Random', s=size, marker="x")
plt.xlim(-0.05,1.05)
# plt.xlim(0.4, 1.01)
# plt.xticks(np.arange(0.4, 1.01, 0.2))
plt.ylim(-0.05,1.05)
plt.xlabel('Recall', fontsize=textsize)
plt.ylabel('Precision', fontsize=textsize)
plt.xticks(fontsize=textsize)
plt.yticks(fontsize=textsize)
plt.legend(bbox_to_anchor=(1.05, 0.4), loc='upper left', fontsize=14, edgecolor='none')
# cbar = plt.colorbar(shrink=0.4, orientation="vertical", drawedges=False, anchor=(1.05, 0.80))
# cbar.ax.tick_params(labelsize=14)
# cbar.ax.set_ylabel('Posterior threshold', fontsize=14, rotation=90)
plt.tight_layout()
plt.savefig(outfile)
plt.close()