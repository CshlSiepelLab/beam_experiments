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
# character_matrix_tsv = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/variable_migration_and_mutation_rates_8_19_24_data_from_8_19_24/raw_data/mig4_mut001_231/mig4_mut001_231_indel_character_matrix.tsv"

outfile = character_matrix_tsv.replace(".tsv", ".pdf")

# Read in the character matrix
df = pd.read_csv(character_matrix_tsv, sep='\t', index_col=0)

# Define color palette
values = np.array([x for sublist in df.values.tolist() for x in sublist])
uniq_values = np.unique(values)
colors = sns.color_palette("hsv", len(uniq_values))

# Map colors to values
color_map = {0: 'white', -1: 'white'}
for i, val in enumerate(uniq_values):
    if val not in color_map:
        color_map[val] = colors[i]

# Plot the character matrix
fig, ax = plt.subplots(figsize=(10, 10))
sns.heatmap(df.values, ax=ax, cmap=list(color_map.values()), cbar=False, xticklabels=False, yticklabels=False, linecolor='black', linewidths=0.5)

# Set the y-axis label and title
ax.set_ylabel("Cells", fontsize=22)
ax.set_title("Barcode mutation sites", fontsize=22)

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(0.5)

# Save the plot to a PDF file
# plt.show()
plt.savefig(outfile)
plt.close()