
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

# Load the data from the CSV file
data = pd.read_csv(infile)

# sort by mmus and then by cp in increasing order
data = data.sort_values(by=["mmus", "cp"])

# make a boxplot of the number of tips per asv for each mmus grouping on the x axis categories
fs = 24
plt.figure(figsize=(14, 8))
strip_plot = sns.stripplot(
    x="mmus", y="information", data=data, color="grey", size=7.5, jitter=True, alpha=0.5
)
box_plot = sns.boxplot(
    x="mmus",
    y="information",
    data=data,
    palette=DEFAULT_COLORS,
    showcaps=True,
    boxprops=dict(facecolor="none", edgecolor="black"),
    whiskerprops=dict(color="black"),
    fliersize=0,
)
unique_mice = data["mmus"].unique()
plt.xlabel("", fontsize=fs)
plt.ylabel("Information", fontsize=fs)
plt.xticks(rotation=45, fontsize=fs)
plt.yticks(fontsize=fs)
# plt.ylim(0, 400)
plt.tight_layout()
plt.savefig(infile.replace(".csv", "_boxplot.pdf"), dpi=1000)
plt.close()
