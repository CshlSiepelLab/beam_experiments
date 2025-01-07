#!/usr/bin/env python3

import pickle as pkl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

file_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/test_metastasis_times.pkl"
origin_time = 250   # given in days
origin_tissue = "P"
min_prob_threshold = 0.5
outfile = "./test_metastasis_times.pdf"

with open(file_path, 'rb') as file:
    met_times = pkl.load(file)

possible_migrations = set()
migration_counts = defaultdict(lambda: np.zeros(origin_time + 1))

prob = 1/len(met_times)

for graph in met_times.values():
    for migration, time in graph.items():
        possible_migrations.add(migration)
        start_range = round(time[0])
        end_range = round(time[1])
        migration_counts[migration][start_range:end_range+1] += prob

df = pd.DataFrame(migration_counts, index=np.arange(0, origin_time + 1)).T

# Split the migration strings into source and target tissues
df.index = pd.MultiIndex.from_tuples([tuple([migration.split('_')[0], "_".join(migration.split('_')[1:])]) for migration in df.index], names=['source', 'target'])

# remove rows without any time over the min probability threshold, indicating no high probability migrations from that source
df = df.loc[df.max(axis=1) >= min_prob_threshold]
remaining_sources = sorted(set(df.index.get_level_values('source').unique()) - {origin_tissue})
remaining_sources.insert(0, origin_tissue)

# Get unique tissues
target_tissues = sorted(set(df.index.get_level_values('target')))
target_tissues_no_num = set([tissue.split('_')[0] for tissue in target_tissues])
tissues = sorted(set(df.index.get_level_values('source')).union(target_tissues_no_num) - {origin_tissue})
tissues.insert(0, origin_tissue)

# Create a grid of subplots with one row per source tissue
fig, axes = plt.subplots(len(remaining_sources), 1, figsize=(15,3 * len(remaining_sources)), sharex=True, sharey=True)
fs = 22

# Create a color palette for the tissues
palette = sns.color_palette("colorblind", len(target_tissues))
target_tissues_reformatted = [tissue.split('_')[0] if "_1" in tissue else tissue for tissue in target_tissues]
tissue_colors = {tissue: palette[i] for i, tissue in enumerate(target_tissues_reformatted)}

for i, source in enumerate(remaining_sources):
    ax = axes[i]
    for target in target_tissues:
        if (source, target) in df.index:
            if "_1" in target:
                target_reformatted = target.split('_')[0]
            else:
                target_reformatted = target
            sns.lineplot(x=df.columns, y=df.loc[(source, target)], ax=ax, color=tissue_colors[target_reformatted])
            ax.fill_between(df.columns, df.loc[(source, target)], alpha=0.3, color=tissue_colors[target_reformatted])
    ax.set_ylabel(source, fontsize=fs)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    if i == len(remaining_sources) - 1:
        ax.set_xlabel('Time (days)', fontsize=fs)

# Create a single legend for all axes
handles = [plt.Line2D([0], [0], color=color, lw=2) for tissue, color in tissue_colors.items()]
labels = list(tissue_colors.keys())
fig.legend(handles, labels, bbox_to_anchor=(1.05, 0.75), title='Target Tissue', frameon=False, fontsize=fs, title_fontsize=fs)

fig.text(0.001, 0.5, 'Source tissue', va='center', ha='center', rotation='vertical', fontsize=fs)

plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust the layout to make room for the legend
plt.savefig(outfile, bbox_inches='tight')
plt.close()
