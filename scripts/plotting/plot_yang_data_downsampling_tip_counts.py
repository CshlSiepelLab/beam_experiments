#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

infile = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/downsampling_yang_2022_real_data_8_22_24/downsampling_tip_counts.csv"
outfile = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/downsampling_yang_2022_real_data_8_22_24/downsampling_tip_counts.pdf"

max_tips=150

# Load the data from the CSV file
data = pd.read_csv(infile)

# For yang data, subtract 1 from the sample number to get the correct sample number without the header
data['Count'] = data['Count'].apply(lambda x: int(x) - 1)

# add a small amount of jitter to the count values to avoid overlapping patches
data['Count'] = data['Count'] + np.random.uniform(0.001, 0.01, size=len(data))

# Create the bar plot
plt.figure(figsize=(14, 8))
bar_plot = sns.barplot(x='Sample', y='Count', hue='Threshold', data=data, palette=DEFAULT_COLORS)

# Add a horizontal line at y=250
plt.axhline(max_tips, color='red', linestyle='--')

# Set the labels and title
plt.xlabel('Sample', fontsize=22)
plt.ylabel('Tip count', fontsize=22)

# Rotate the x-axis labels for better readability
plt.xticks(rotation=90, fontsize=22)
plt.yticks(fontsize=22)

# Add a legend for the thresholds
# plt.legend(title="Downsampling\ndistance\nthreshold", fontsize=22, title_fontsize=22, bbox_to_anchor=(1.05, 1), loc='upper left')

# Add outline to the first bar in each x-axis group that is <= 250
grouped = data.groupby('Sample')
for name, group in grouped:
    first_below_threshold = group[(group['Count'] <= max_tips) & (group['Threshold'] >= 0.1)].head(1)
    print(first_below_threshold)
    if first_below_threshold.empty:
        first_below_threshold = group.sort_values(by='Threshold', ascending=False).head(1)
    index = first_below_threshold.index[0]

    for i, patch in enumerate(bar_plot.patches):
        if patch.get_height() == first_below_threshold['Count'].values[0]:
            index = i
            break

    bar_plot.patches[index].set_edgecolor('black')
    bar_plot.patches[index].set_linewidth(2)

# Add a custom legend entry for the black box
handles, labels = bar_plot.get_legend_handles_labels()
handles.append(plt.Line2D([0], [0], color='black', linewidth=2))
labels.append('Selected threshold')

# Update the legend with the new entry
plt.legend(handles=handles, labels=labels, title="Downsampling\ndistance\nthreshold", fontsize=22, title_fontsize=22, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

# Show the plot
plt.tight_layout()
plt.savefig(outfile, dpi=2000)
plt.close()