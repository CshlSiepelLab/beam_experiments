#!/usr/bin/env python

import sys
import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# inputs
character_matrix_tsv = sys.argv[1]

# # testing
# character_matrix_tsv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1544/CP01/downsampled_char_matrix.tsv"

outfile = character_matrix_tsv.replace(".tsv", ".pdf")

# Read in the character matrix
df = pd.read_csv(character_matrix_tsv, sep='\t', index_col=0)

# Define color palette
uniq_values = np.unique(np.array([x for sublist in df.values.tolist() for x in sublist]))

# Map colors to values
color_map = {'0': '#FFFFFF', '-1': '#808080'}
uniq_values = uniq_values[~np.isin(uniq_values, [0, -1])]
for val in uniq_values:
    color = '#%06X' % np.random.randint(0, 0xFFFFFF)
    while color in color_map.values():
        color = '#%06X' % np.random.randint(0, 0xFFFFFF)
    color_map[str(val)] = color

# Plot the character matrix
fig, ax = plt.subplots(figsize=(10, 10))

sns.heatmap(df.values, ax=ax, cbar=False, xticklabels=False, yticklabels=False, linecolor='black', linewidths=0.5)

# Iterate over the cells and change the color based on the value
for i in range(df.shape[0]):
    for j in range(df.shape[1]):
        ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, edgecolor='black', facecolor=color_map.get(str(df.values[i, j]), '#FFFFFF')))

# Set the y-axis label and title
ax.set_ylabel("Cells", fontsize=32)
ax.set_title("Barcode mutation sites", fontsize=32)

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(0.5)

# Save the plot to a PDF file
# plt.show()
plt.savefig(outfile)
plt.close()