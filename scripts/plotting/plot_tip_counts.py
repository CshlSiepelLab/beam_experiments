#!/usr/bin/env python3

import sys
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

infile = sys.argv[1]
# infile = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/asv_counts_per_cp.csv"

max_tips = int(sys.argv[2])
# max_tips = 150

# Load the data from the CSV file
data = pd.read_csv(infile)

# sort by mouse and then by cp in increasing order
data = data.sort_values(by=["mouse", "cp"])

# remove any CP00 samples
data = data[data["cp"] != "CP0"]
data = data[data["cp"] != "CP00"]
data = data[data["cp"] != "CP000"]

# Create the bar plot
fs = 22
plt.figure(figsize=(14, 8))
bar_plot = sns.barplot(
    x="mouse",
    y="num_asvs",
    hue="cp",
    width=0.95,
    data=data,
    ci=None,
    palette=DEFAULT_COLORS,
)
bar_plot.legend_.remove()
# plt.axhline(max_tips, color='red', linestyle='--')
plt.xlabel("Sample clonal populations", fontsize=22)
plt.ylabel("Tip count", fontsize=22)
plt.xticks(rotation=90, fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
plt.savefig(infile.replace(".csv", "_barplot.pdf"), dpi=1000)
plt.close()

# make a boxplot of the number of tips per asv for each mouse grouping on the x axis categories
fs = 22
plt.figure(figsize=(14, 8))
strip_plot = sns.stripplot(
    x="mouse", y="num_asvs", data=data, color="grey", size=7.5, jitter=True, alpha=0.5
)
box_plot = sns.boxplot(
    x="mouse",
    y="num_asvs",
    data=data,
    palette=DEFAULT_COLORS,
    showcaps=True,
    boxprops=dict(facecolor="none", edgecolor="black"),
    whiskerprops=dict(color="black"),
    fliersize=0,
)
plt.axhline(max_tips, color="red", linestyle="--")
unique_mice = data["mouse"].unique()
for i, mouse in enumerate(unique_mice):
    mouse_data = data[data["mouse"] == mouse]
    for j in range(mouse_data.shape[0]):
        if mouse_data.iloc[j]["num_asvs"] >= max_tips:
            plt.text(
                i,
                mouse_data.iloc[j]["num_asvs"] + 5,
                mouse_data.iloc[j]["cp"],
                color="black",
                ha="center",
                fontsize=fs,
            )
plt.xlabel("", fontsize=fs)
plt.ylabel("Tip count per clonal population", fontsize=fs)
plt.xticks(rotation=45, fontsize=fs)
plt.yticks(fontsize=fs)
# plt.ylim(0, 200)
plt.tight_layout()
plt.savefig(infile.replace(".csv", "_boxplot.pdf"), dpi=1000)
plt.close()
