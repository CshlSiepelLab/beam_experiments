#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

infile='/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_10_31_24/downsampled_data/tip_counts.csv'
outfile='/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_10_31_24/downsampled_data/tip_counts.pdf'

# Read the CSV file into a DataFrame
df = pd.read_csv(infile)

# Convert tip_count to percentages
df['tip_count'] = df['tip_count'] / 100

# Create a bar plot with points
plt.figure(figsize=(10, 6))
sns.boxplot(x='downsampling_threshold', y='tip_count', data=df, palette=DEFAULT_COLORS, boxprops=dict(alpha=0.3), showfliers=False)
sns.stripplot(x='downsampling_threshold', y='tip_count', data=df, hue='downsampling_threshold', palette=DEFAULT_COLORS, jitter=True, size=5, legend=False)

# Set plot labels and title
plt.xlabel('Downsampling distance threshold', fontsize=26)
plt.ylabel('% original tips retained', fontsize=26)
plt.xticks(fontsize=26)
plt.yticks(fontsize=26)
plt.tight_layout()

plt.savefig(outfile)
plt.close()

