#!/usr/bin/env python3

import pandas as pd
import ast
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


def calculate_f1_score(precision, recall):
    if precision + recall == 0:
        return 0
    return 2 * (precision * recall) / (precision + recall)


def get_number_of_migrations_from_graph_dict(graph_dict):
    graph_dict = ast.literal_eval(graph_dict)
    return sum(graph_dict.values())


# Load the data
csv_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_tree_test_100_tips_10_31_24/precision_recall/downsample_all_stats_formatted.csv"
data = pd.read_csv(csv_file)

tip_count_csv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_tree_test_100_tips_10_31_24/precision_recall/tip_counts.csv"
tip_count_data = pd.read_csv(tip_count_csv)

# Merge the data with the tip count data
data = pd.merge(
    data,
    tip_count_data,
    left_on=["sim", "downsample_threshold"],
    right_on=["sim_name", "downsampling_threshold"],
)

# Replace the f1_score calculation with the custom function
data["f1_score"] = data.apply(
    lambda row: calculate_f1_score(row["precision"], row["recall"]), axis=1
)

# get the migration counts
data["migration_count"] = data.apply(
    lambda row: get_number_of_migrations_from_graph_dict(row["thresh_counts"]), axis=1
)

# Plotting
fontsize = 24
# plt.figure(figsize=(10, 6))
# for sim, group in data.groupby('sim'):
#     plt.plot(group['migration_count'], group['f1_score'], marker='o')

# plt.xlabel('Downsample Threshold')
# plt.ylabel('F1 Score')
# plt.savefig(csv_file.replace(".csv", "_trace_f1.pdf"))
# plt.close()

# filter data to only include the <= a specified downsample threshold
data = data[data["downsample_threshold"] <= 1.0]

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
unique_thresholds = sorted(data["downsample_threshold"].unique())
color_map = {
    threshold: DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
    for i, threshold in enumerate(unique_thresholds)
}


# plot tip count by migration count
plt.figure(figsize=(10, 6))
for sim, group in data.groupby("sim"):
    plt.plot(
        group["tip_count"],
        group["migration_count"],
        color="lightgrey",
        linestyle="-",
        linewidth=1,
        zorder=1,
    )
    scatter = plt.scatter(
        group["tip_count"],
        group["migration_count"],
        c=group["downsample_threshold"].map(color_map),
        marker="o",
        zorder=2,
    )

# Calculate Pearson correlation coefficient
pearson_corr, _ = pearsonr(data["tip_count"], data["migration_count"])

# Draw the linear regression line
m, b = np.polyfit(data["tip_count"], data["migration_count"], 1)
plt.plot(
    data["tip_count"],
    m * data["tip_count"] + b,
    color="black",
    linestyle="--",
    linewidth=2,
)

# Create a legend for the colors
handles = [
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10)
    for threshold, color in color_map.items()
]
labels = [f"{threshold}" for threshold in color_map.keys()]
plt.legend(
    handles,
    labels,
    title="Downsampling\ndistance\nthreshold",
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
    frameon=False,
    fontsize=16,
    title_fontsize=16,
)

# Add the Pearson correlation coefficient to the title
plt.text(
    0.65,
    0.95,
    f"Pearson r={pearson_corr:.2f}",
    fontsize=16,
    transform=plt.gca().transAxes,
    verticalalignment="top",
)
plt.text(
    0.05,
    0.95,
    f"y={m:.2f}x + {b:.2f}",
    fontsize=16,
    transform=plt.gca().transAxes,
    verticalalignment="top",
)
plt.ylim(0, 70)
plt.xlabel("Leaf count", fontsize=fontsize)
plt.ylabel("Migration count", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.tight_layout()
plt.savefig(csv_file.replace(".csv", "_tip_count_by_migration_count.pdf"))
plt.close()
