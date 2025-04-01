#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import re
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt

# input_file = sys.argv[1]
# primary_tissue = sys.argv[2]

input_file = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1467/CP01/combined.log"
primary_tissue = "PRL"

outfile = input_file.replace(".log", "_tissue_rate_matrix_hypothesis_setups.pdf")

df = pd.read_csv(input_file, sep="\t", comment="#")

tissue_rate_col_names = [
    name for name in df.columns if name.startswith("geoSubstModelLogger")
]
tissues = list(
    set(
        [
            tissue
            for name in tissue_rate_col_names
            for tissue in name.replace("geoSubstModelLogger.relGeoRate_", "").split("_")
        ]
    )
)
tissues = [primary_tissue] + sorted([tis for tis in tissues if tis != primary_tissue])

# Create a matrix for the rates with labels
num_tissues = len(tissues)
rate_matrix_labels = np.empty((num_tissues, num_tissues), dtype=object)
rate_matrix_primary_zero = np.empty((num_tissues, num_tissues), dtype=object)

for i, source in enumerate(tissues):
    for j, recipient in enumerate(tissues):
        if source == recipient:
            rate_matrix_labels[i, j] = -1
            rate_matrix_primary_zero[i, j] = -1
        else:
            if recipient == primary_tissue:
                rate_matrix_primary_zero[i, j] = 0
                rate_matrix_labels[i, j] = 3  # set specific to 4x4 matrix
            elif source == primary_tissue:
                rate_matrix_primary_zero[i, j] = 1
                rate_matrix_labels[i, j] = 1
            else:
                rate_matrix_primary_zero[i, j] = 2
                rate_matrix_labels[i, j] = 2

rate_matrix_labels = rate_matrix_labels.astype(int)
rate_matrix_primary_zero = rate_matrix_primary_zero.astype(int)

# Plot the rate matrices
fs = 36
fig, axes = plt.subplots(1, 2, figsize=(20, 8))

colors = ["grey", "white", "blue", "orange"]
# colors = ['grey', 'white', "#6fa8dc", "#e69138"]
cmap_null = ListedColormap(colors)

# colors = ['grey', 'white', "#6fa8dc", "#e69138", "#c27ba0"]
colors = ["grey", "white", "blue", "orange", "green"]
cmap_alt = ListedColormap(colors)


sns.heatmap(
    rate_matrix_primary_zero,
    xticklabels=tissues,
    yticklabels=tissues,
    annot=True,
    fmt="",
    cmap=cmap_null,
    cbar=False,
    ax=axes[0],
    annot_kws={"size": fs},
)
axes[0].set_title("H$_{null}$", fontsize=fs)
axes[0].set_xlabel("Recipient", fontsize=fs)
axes[0].set_ylabel("Source", fontsize=fs)
axes[0].tick_params(labelsize=fs)

sns.heatmap(
    rate_matrix_labels,
    xticklabels=tissues,
    yticklabels=tissues,
    annot=True,
    fmt="",
    cmap=cmap_alt,
    cbar=False,
    ax=axes[1],
    annot_kws={"size": fs},
)
axes[1].set_title("H$_{alt}$", fontsize=fs)
axes[1].set_xlabel("Recipient", fontsize=fs)
axes[1].set_ylabel("Source", fontsize=fs)
axes[1].tick_params(labelsize=fs)

plt.tight_layout()
plt.savefig(outfile)
plt.close()
