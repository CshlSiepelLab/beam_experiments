#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

DEFAULT_COLORS = [
    "#6aa84f",
    "#be5742e1",
    "#6fa8dc",
    "#e69138",
    "#9e9e9e",
    "#c27ba0",
    "brown",
    "black",
    "darkgreen",
    "purple",
    "blue",
] * 3

infile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsample_test_yang_2022_real_data_10_1_24/downsampling_tip_counts.csv"
outfile = infile.replace(".csv", ".pdf")

max_tips = 250

# Load the data from the CSV file
data = pd.read_csv(infile)

# For yang data, subtract 1 from the sample number to get the correct sample number without the header
data["Count"] = data["Count"].apply(lambda x: int(x) - 1)

# add a small amount of jitter to the count values to avoid overlapping patches
data["Count"] = data["Count"] + np.random.uniform(0.001, 0.01, size=len(data))

# sort by genotype
data["Genotype"] = data["Sample"].apply(lambda x: x.split("_")[1])
data = data.sort_values(by="Genotype")

# Create the bar plot
plt.figure(figsize=(14, 8))
bar_plot = sns.barplot(
    x="Sample", y="Count", hue="Threshold", data=data, palette=DEFAULT_COLORS
)

# Add a horizontal line at y=250
# plt.axhline(max_tips, color='red', linestyle='--')

# Set the labels and title
plt.xlabel("Sample", fontsize=22)
plt.ylabel("Tip count", fontsize=22)

# Rotate the x-axis labels for better readability
plt.xticks(rotation=90, fontsize=22)
plt.yticks(np.arange(0, data["Count"].max() + 200, 200), fontsize=22)

# Add a legend for the thresholds
plt.legend(
    title="Downsampling\ndistance\nthreshold",
    fontsize=22,
    title_fontsize=22,
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
    frameon=False,
)

#### OUTLINE BELOW NOT WORKING PROPERLY
# # Add outline to the first bar in each x-axis group that is <= max tips and >= 0.1 threshold
# grouped = data.groupby('Sample')

# for name, group in grouped:
#     first_below_threshold = group[(group['Count'] <= max_tips) & (group['Threshold'] >= 0.1)].head(1)
#     if first_below_threshold.empty:
#         first_below_threshold = group.sort_values(by='Threshold', ascending=False).head(1)

#     # Find the height of the selected bar
#     selected_height = first_below_threshold['Count'].values[0]

#     # Find the correct patch based on height
#     for patch in bar_plot.patches:
#         if patch.get_height() == selected_height:
#             patch.set_edgecolor('black')
#             patch.set_linewidth(2)
#             break

# # Add a custom legend entry for the black box
# handles, labels = bar_plot.get_legend_handles_labels()
# handles.append(plt.Line2D([0], [0], color='black', linewidth=2))
# labels.append('Selected threshold')

# Update the legend with the new entry
# plt.legend(handles=handles, labels=labels, title="Downsampling\ndistance\nthreshold", fontsize=22, title_fontsize=22, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

# Show the plot
plt.tight_layout()
plt.savefig(outfile, dpi=2000)
plt.close()
