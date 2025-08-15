
import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ete3 import Tree
import dendropy
from copy import deepcopy
import ast
import pickle


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


def posterior_threshold_metrics(posterior_prob_graph, true_counts, i, t=0.50):
    max_threshold = max(map(float, posterior_prob_graph.values()))
    thresholds = [j for j in np.arange(0, max_threshold, 0.01)]
    rows = []
    for thresh in thresholds:
        thresh_counts = {
            key: value for key, value in posterior_prob_graph.items() if value > thresh
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

    # check that the max for the data is not below the threshold for the consensus graph
    if max_threshold < t:
        t = np.floor(max_threshold * 100) / 100
    threshold_df = thresh_prec_rec[np.isclose(thresh_prec_rec["Threshold"], t)]
    post_prob_precision = threshold_df["precision"].values[0]
    post_prob_recall = threshold_df["recall"].values[0]
    post_prob_f1 = 2 * (
        (post_prob_precision * post_prob_recall)
        / (post_prob_precision + post_prob_recall)
    )

    return thresh_prec_rec, rows, post_prob_f1, post_prob_recall, post_prob_precision


# user inputs
dirs = (sys.argv[1]).split(",")
primary_tissue = sys.argv[2]
outdir = sys.argv[3]
plot_independent_plots = False

# make a file to record performance statistics for all sim datasets
outfile_metrics = f"{outdir}/metrics.csv"
with open(outfile_metrics, "w") as file:
    header = "sim,Random_f1,Random_recall,Random_precision,Consensus_f1,Consensus_recall,Consensus_precision,Parsimony_f1,Parsimony_recall,Parsimony_precision,MACHINA_f1,MACHINA_recall,MACHINA_precision,MACH2_f1,MACH2_recall,MACH2_precision,Metient_f1,Metient_recall,Metient_precision,BEAM_f1,BEAM_recall,BEAM_precision\n"
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
mach2_all_thresh_rows = []
all_thresh_rows = []
i = 0

for true_tree_file in dirs:

    # makes sure true_tree_file is not an empty string
    if not true_tree_file:
        continue

    # get working dir from outdir based on snakemake input
    dir = os.path.dirname(os.path.dirname(os.path.dirname(true_tree_file)))
    sim = os.path.basename(os.path.dirname(true_tree_file))

    # get other files to compare
    machina_file = f"{dir}/machina/{sim}/machina_tree_all_tissue_labels.nwk"
    mach2_file = f"{dir}/mach2/{sim}/consensus_graph.txt"
    metient_file = f"{dir}/metient/{sim}/{sim}_{primary_tissue}_migration_graphs.txt"
    
    # metient_file=f"/grid/siepel/home/staklins/stored_results/beam/latest_results/snakemake_performance_uniform_50cells_50sites_data_7_24_24/metient_calibrate_80_ideal_sims_8_15_25/calibrate/{sim}_migration_graphs.txt"
    
    beast_posterior_file = f"{dir}/beam_gtr/{sim}/posterior_prob_graph.csv"
    consensus_file = (
        f"{dir}/random_consensus_parsimony_tissue_inference/{sim}/consensus_tissues.nwk"
    )
    random_file = (
        f"{dir}/random_consensus_parsimony_tissue_inference/{sim}/random_tissues.nwk"
    )
    parsimony_file = (
        f"{dir}/random_consensus_parsimony_tissue_inference/{sim}/parsimony_tissues.nwk"
    )
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

    # output true graph to a file
    true_graph_dir = f"{outdir}/true_graphs"
    if not os.path.exists(true_graph_dir):
        os.makedirs(true_graph_dir)
    with open(f"{true_graph_dir}/{sim}_true_graph.csv", "w") as file:
        file.write(f"source_target,num_edges\n")
        for key, value in true_counts.items():
            file.write(f"{key},{value}\n")

    # process machina result to get precision and recall
    machina_f1 = float("nan")
    machina_recall = float("nan")
    machina_precision = float("nan")
    if os.path.exists(machina_file):
        machina_counts = process_tree(machina_file)
        machina_f1, machina_recall, machina_precision = calculate_metrics(
            true_counts, machina_counts
        )
        machina_precisions[i] = machina_precision
        machina_recalls[i] = machina_recall

    # process metient result to get precision and recall
    metient_f1 = float("nan")
    metient_recall = float("nan")
    metient_precision = float("nan")
    if os.path.exists(metient_file):
        with open(metient_file, "r") as file:
            lines = file.readlines()
        top_loss_metient = lines[1].strip().split("\t")
        metient_counts_input = ast.literal_eval(top_loss_metient[1])
        metient_counts = {}
        for outer_key, inner_dict in metient_counts_input.items():
            for inner_key, value in inner_dict.items():
                if value != 0.0:
                    metient_counts[f"{outer_key}_{inner_key}"] = value
        metient_f1, metient_recall, metient_precision = calculate_metrics(
            true_counts, metient_counts
        )
        metient_precisions[i] = metient_precision
        metient_recalls[i] = metient_recall

    # process random result to get precision and recall
    random_f1 = float("nan")
    random_recall = float("nan")
    random_precision = float("nan")
    if os.path.exists(random_file):
        random_counts = process_tree(random_file)
        random_f1, random_recall, random_precision = calculate_metrics(
            true_counts, random_counts
        )
        random_precisions[i] = random_precision
        random_recalls[i] = random_recall

    # process consensus result to get precision and recall
    consensus_f1 = float("nan")
    consensus_recall = float("nan")
    consensus_precision = float("nan")
    if os.path.exists(consensus_file):
        consensus_counts = process_tree(consensus_file)
        consensus_f1, consensus_recall, consensus_precision = calculate_metrics(
            true_counts, consensus_counts
        )
        consensus_precisions[i] = consensus_precision
        consensus_recalls[i] = consensus_recall

    # process greedy fitch parsimony result to get precision and recall
    parsimony_f1 = float("nan")
    parsimony_recall = float("nan")
    parsimony_precision = float("nan")
    if os.path.exists(parsimony_file):
        parsimony_counts = process_tree(parsimony_file)
        parsimony_f1, parsimony_recall, parsimony_precision = calculate_metrics(
            true_counts, parsimony_counts
        )
        parsimony_precisions[i] = parsimony_precision
        parsimony_recalls[i] = parsimony_recall

    # process mach2 result
    mach2_f1 = float("nan")
    mach2_recall = float("nan")
    mach2_precision = float("nan")
    if os.path.exists(mach2_file):
        mach2_prob_graph = {}
        with open(mach2_file, "r") as file:
            for line in file.readlines():
                line = line.strip().split(",")
                mach2_prob_graph[line[0]] = float(line[1])
        mach2_thresh_prec_rec, mach2_rows, mach2_f1, mach2_recall, mach2_precision = (
            posterior_threshold_metrics(mach2_prob_graph, true_counts, sim)
        )
        mach2_all_thresh_rows.extend(mach2_rows)

    # process beast posterior
    post_prob_f1 = float("nan")
    post_prob_recall = float("nan")
    post_prob_precision = float("nan")
    if os.path.exists(beast_posterior_file):
        posterior_prob_graph = {}
        with open(beast_posterior_file, "r") as file:
            for line in file.readlines():
                line = line.strip().split(",")
                posterior_prob_graph[line[0]] = float(line[1])
        thresh_prec_rec, rows, post_prob_f1, post_prob_recall, post_prob_precision = (
            posterior_threshold_metrics(posterior_prob_graph, true_counts, sim)
        )
        all_thresh_rows.extend(rows)

    # write metrics used for the plot to a file
    with open(outfile_metrics, "a") as file:
        data = f"{sim},{random_f1},{random_recall},{random_precision},{consensus_f1},{consensus_recall},{consensus_precision},{parsimony_f1},{parsimony_recall},{parsimony_precision},{machina_f1},{machina_recall},{machina_precision},{mach2_f1},{mach2_recall},{mach2_precision},{metient_f1},{metient_recall},{metient_precision},{post_prob_f1},{post_prob_recall},{post_prob_precision}\n"
        file.write(data)

    if plot_independent_plots:
        size = 75
        textsize = 18
        plt.figure()
        if os.path.exists(beast_posterior_file):
            # plt.scatter(thresh_prec_rec['recall'], thresh_prec_rec['precision'], c=thresh_prec_rec['Threshold'], cmap='viridis', s=25, marker='x')
            plt.plot(
                thresh_prec_rec["recall"],
                thresh_prec_rec["precision"],
                color="red",
                label="BEAM",
            )
        if os.path.exists(mach2_file):
            # plt.scatter(mach2_thresh_prec_rec['recall'], mach2_thresh_prec_rec['precision'], c=mach2_thresh_prec_rec['Threshold'], cmap='viridis', s=25, marker='x')
            plt.plot(
                mach2_thresh_prec_rec["recall"],
                mach2_thresh_prec_rec["precision"],
                color="navy",
                label="Mach2",
            )
        if os.path.exists(machina_file):
            plt.scatter(
                machina_recall,
                machina_precision,
                color="gold",
                label="Machina",
                s=size,
                marker="x",
            )
        if os.path.exists(metient_file):
            plt.scatter(
                metient_recall,
                metient_precision,
                color="green",
                label="Metient",
                s=size,
                marker="x",
            )
        if os.path.exists(consensus_file):
            plt.scatter(
                consensus_recall,
                consensus_precision,
                color="blue",
                label="Consensus",
                s=size,
                marker="x",
            )
        if os.path.exists(random_file):
            plt.scatter(
                random_recall,
                random_precision,
                color="black",
                label="Random",
                s=size,
                marker="x",
            )
        if os.path.exists(parsimony_file):
            plt.scatter(
                parsimony_recall,
                parsimony_precision,
                color="Purple",
                label="Parsimony",
                s=size,
                marker="x",
            )
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)
        plt.xlabel("Recall", fontsize=textsize)
        plt.ylabel("Precision", fontsize=textsize)
        plt.xticks(fontsize=textsize)
        plt.yticks(fontsize=textsize)
        plt.legend(
            bbox_to_anchor=(1.05, 0.4), loc="upper left", fontsize=14, edgecolor="none"
        )
        # cbar = plt.colorbar(shrink=0.4, orientation="vertical", drawedges=False, anchor=(1.05, 0.80))
        # cbar.ax.tick_params(labelsize=14)
        # cbar.ax.set_ylabel('Posterior threshold', fontsize=14, rotation=90)
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

    i += 1

mach2_all_thresh_df = pd.DataFrame(mach2_all_thresh_rows)
all_thresh_df = pd.DataFrame(all_thresh_rows)

# save intermediate variables to a file for later re-plotting
with open(f"{outdir}/precision_recall_vars.pkl", "wb") as file:
    pickle.dump(
        [
            machina_precisions,
            machina_recalls,
            metient_precisions,
            metient_recalls,
            random_precisions,
            random_recalls,
            consensus_precisions,
            consensus_recalls,
            parsimony_precisions,
            parsimony_recalls,
            mach2_all_thresh_df,
            all_thresh_df,
        ],
        file,
    )

# # optionally open from pickle file and avoid recalculations above
# with open(f"{outdir}/precision_recall_vars.pkl", "rb") as file:
#     machina_precisions, machina_recalls, metient_precisions, metient_recalls, random_precisions, random_recalls, consensus_precisions, consensus_recalls, parsimony_precisions, parsimony_recalls, mach2_all_thresh_df, all_thresh_df = pickle.load(file)

# make an overall averaged precision/recall curve
outfile = f"{outdir}/precision_recall.pdf"
avg_machina_precision = float("nan")
avg_machina_recall = float("nan")
avg_metient_precision = float("nan")
avg_metient_recall = float("nan")
avg_random_precision = float("nan")
avg_random_recall = float("nan")
avg_consensus_precision = float("nan")
avg_consensus_recall = float("nan")
avg_parsimony_precision = float("nan")
avg_parsimony_recall = float("nan")


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

if not all_thresh_df.empty:
    avg_df = (
        all_thresh_df.groupby("Threshold")[["precision", "recall"]].mean().reset_index()
    )
    all_thresh_df.to_csv(f"{outdir}/beam_all_threshold_stats.csv", index=False)

if not mach2_all_thresh_df.empty:
    avg_mach2_df = (
        mach2_all_thresh_df.groupby("Threshold")[["precision", "recall"]]
        .mean()
        .reset_index()
    )
    mach2_all_thresh_df.to_csv(f"{outdir}/mach2_all_threshold_stats.csv", index=False)

size = 100
textsize = 20
plt.figure()
if not avg_df.empty:
    # plt.scatter(avg_df['recall'], avg_df['precision'], c=avg_df['Threshold'], cmap='viridis', s=25, marker='x')
    plt.plot(avg_df["recall"], avg_df["precision"], color="red", label="BEAM")
if not avg_mach2_df.empty:
    # plt.scatter(avg_mach2_df['recall'], avg_mach2_df['precision'], c=avg_mach2_df['Threshold'], cmap='viridis', s=25, marker='x')
    plt.plot(
        avg_mach2_df["recall"], avg_mach2_df["precision"], color="navy", label="MACH2"
    )
if not np.isnan(avg_machina_recall) and not np.isnan(avg_machina_precision):
    plt.scatter(
        avg_machina_recall,
        avg_machina_precision,
        color="gold",
        label="MACHINA",
        s=size,
        marker="x",
    )
if not np.isnan(avg_metient_recall) and not np.isnan(avg_metient_precision):
    plt.scatter(
        avg_metient_recall,
        avg_metient_precision,
        color="green",
        label="Metient",
        s=size,
        marker="x",
    )
if not np.isnan(avg_consensus_recall) and not np.isnan(avg_consensus_precision):
    plt.scatter(
        avg_consensus_recall,
        avg_consensus_precision,
        color="blue",
        label="Consensus",
        s=size,
        marker="x",
    )
if not np.isnan(avg_random_recall) and not np.isnan(avg_random_precision):
    plt.scatter(
        avg_random_recall,
        avg_random_precision,
        color="black",
        label="Random",
        s=size,
        marker="x",
    )
if not np.isnan(avg_parsimony_recall) and not np.isnan(avg_parsimony_precision):
    plt.scatter(
        avg_parsimony_recall,
        avg_parsimony_precision,
        color="purple",
        label="Parsimony",
        s=size,
        marker="x",
    )
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.xlabel("Recall", fontsize=textsize)
plt.ylabel("Precision", fontsize=textsize)
plt.xticks(fontsize=textsize)
plt.yticks(fontsize=textsize)
plt.legend(bbox_to_anchor=(1.05, 0.5), loc="upper left", fontsize=14, edgecolor="none")
plt.tight_layout()
plt.savefig(outfile)
plt.close()
