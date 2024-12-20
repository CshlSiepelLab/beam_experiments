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
import glob
import ast
import pickle

# meant to be used after the initial precision/recall processing across all datasets, now stratifying them by topology restriction in the simulation
# outdir = sys.argv[1]

# testing
outdir = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/precision_recall_curve"

# optionally to open from pickle file and avoid recalculations above
with open(f"{outdir}/precision_recall_vars.pkl", "rb") as file:
    # machina_precisions, machina_recalls, metient_precisions, metient_recalls, random_precisions, random_recalls, consensus_precisions, consensus_recalls, parsimony_precisions, parsimony_recalls, all_thresh_df = pickle.load(file)
    machina_precisions, machina_recalls, metient_precisions, metient_recalls, random_precisions, random_recalls, consensus_precisions, consensus_recalls, parsimony_precisions, parsimony_recalls, all_thresh_df, pathfinder_all_thresh_df = pickle.load(file)

all_thresh_df['sim_topology'] = all_thresh_df['sim'].apply(lambda x: x.split('_')[0])

# topologies = all_thresh_df['sim_topology'].unique()
topologies = ['mS', 'pS', 'pM', 'pR']
num_classes = len(topologies)

sim_names = all_thresh_df['sim'].unique()


# Group sims by topology based on index assuming sim_names is the order of values in all other lists
topology_groups = {topology: [] for topology in topologies}
for sim_name in sim_names:
    topology = sim_name.split('_')[0]
    topology_groups[topology].append(sim_name)

# Create 4 axes plots, one for each topology
fig, axes = plt.subplots(1, num_classes, figsize=(5 * num_classes, 5), sharey=False)
axes = axes.flatten()
size = 100
textsize = 22

for i, (topology, sims) in enumerate(topology_groups.items()):
    ax = axes[i]
    topology_df = all_thresh_df[all_thresh_df['sim'].isin(sims)]
    avg_df = topology_df.groupby('Threshold')[['precision', 'recall']].mean().reset_index()

    avg_machina_precision = np.nanmean([machina_precisions[sim_names.tolist().index(sim)] for sim in sims])
    avg_machina_recall = np.nanmean([machina_recalls[sim_names.tolist().index(sim)] for sim in sims])
    avg_metient_precision = np.nanmean([metient_precisions[sim_names.tolist().index(sim)] for sim in sims])
    avg_metient_recall = np.nanmean([metient_recalls[sim_names.tolist().index(sim)] for sim in sims])
    avg_random_precision = np.nanmean([random_precisions[sim_names.tolist().index(sim)] for sim in sims])
    avg_random_recall = np.nanmean([random_recalls[sim_names.tolist().index(sim)] for sim in sims])
    avg_consensus_precision = np.nanmean([consensus_precisions[sim_names.tolist().index(sim)] for sim in sims])
    avg_consensus_recall = np.nanmean([consensus_recalls[sim_names.tolist().index(sim)] for sim in sims])
    avg_parsimony_precision = np.nanmean([parsimony_precisions[sim_names.tolist().index(sim)] for sim in sims])
    avg_parsimony_recall = np.nanmean([parsimony_recalls[sim_names.tolist().index(sim)] for sim in sims])

    if not avg_df.empty:
        ax.plot(avg_df['recall'], avg_df['precision'], color='grey', label='BEAM')
    if not np.isnan(avg_machina_recall) and not np.isnan(avg_machina_precision):
        ax.scatter(avg_machina_recall, avg_machina_precision, color='red', label='MACHINA', s=size, marker="x")
    if not np.isnan(avg_metient_recall) and not np.isnan(avg_metient_precision):
        ax.scatter(avg_metient_recall, avg_metient_precision, color='green', label='Metient', s=size, marker="x")
    if not np.isnan(avg_consensus_recall) and not np.isnan(avg_consensus_precision):
        ax.scatter(avg_consensus_recall, avg_consensus_precision, color='blue', label='Consensus', s=size, marker="x")
    if not np.isnan(avg_random_recall) and not np.isnan(avg_random_precision):
        ax.scatter(avg_random_recall, avg_random_precision, color='black', label='Random', s=size, marker="x")
    if not np.isnan(avg_parsimony_recall) and not np.isnan(avg_parsimony_precision):
        ax.scatter(avg_parsimony_recall, avg_parsimony_precision, color='purple', label='Parsimony', s=size, marker="x")
    
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('Recall', fontsize=textsize)
    ax.set_ylabel('Precision', fontsize=textsize)
    ax.set_title(f'{topology}', fontsize=textsize)
    ax.tick_params(axis='both', which='major', labelsize=textsize)
    if i == len(topologies) - 1:
        ax.legend(bbox_to_anchor=(1.05, 0.5), loc='upper left', fontsize=14, edgecolor='none')

plt.tight_layout()
# plt.subplots_adjust(left=0)
outfile = f"{outdir}/precision_recall_by_seeding_topology.pdf"
plt.savefig(outfile, bbox_inches="tight")
plt.close()
