#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

infile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/asv_counts_per_cp.csv"
outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/asv_counts_per_cp.pdf"

max_tips=200

# Load the data from the CSV file
data = pd.read_csv(infile)

# Create the bar plot
plt.figure(figsize=(14, 8))
bar_plot = sns.barplot(x='mouse', y='num_asvs', hue='cp', width = 0.95, data=data, ci=None, palette=DEFAULT_COLORS)

# turn off legend
bar_plot.legend_.remove()

# Add a horizontal line
# plt.axhline(max_tips, color='red', linestyle='--')

# Set the labels and title
plt.xlabel('Sample clonal populations', fontsize=22)
plt.ylabel('Tip count', fontsize=22)

# Rotate the x-axis labels for better readability
plt.xticks(rotation=90, fontsize=22)
plt.yticks(fontsize=22)

# Show the plot
plt.tight_layout()
plt.savefig(outfile, dpi=1000)
plt.close()