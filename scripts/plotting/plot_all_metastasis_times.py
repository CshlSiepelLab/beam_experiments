
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


file_path = sys.argv[1]
origin_time = float(sys.argv[2])
origin_tissue = sys.argv[3]
outfile = sys.argv[4]


data = pd.read_csv(file_path)

# Replace underscores in source_target for better readability
data["source_target"] = data["source_target"].str.replace("_", "->")

# Group by source tissue and then sort by time within each source tissue
data["source"] = data["source_target"].apply(lambda x: x.split("->")[0])
data["target"] = data["source_target"].apply(lambda x: x.split("->")[1])
data["mean_time"] = data.groupby(["source_target", "threshold"])["time"].transform(
    "mean"
)
data = data.sort_values(by=["mean_time", "source"])

fs = 14

# Create a grid of plots, one for each threshold
g = sns.FacetGrid(
    data=data, col="threshold", col_wrap=3, height=4, sharex=True, sharey=True
)

# Use a boxplot without outliers and overlay data points
g.map(
    sns.boxplot,
    "time",
    "source_target",
    order=data["source_target"].unique(),
    palette="tab10",
    orient="h",
    showfliers=False,
)
g.map(
    sns.stripplot,
    "time",
    "source_target",
    order=data["source_target"].unique(),
    palette="dark:.3",
    orient="h",
    size=3,
    jitter=True,
    alpha=0.6,
)

# Add a title and labels
g.set_axis_labels("Time", "Source -> Target")
g.set(xlim=(0, origin_time))

# Adjust text sizes
g.set_titles("{col_name} threshold", size=fs)
g.set_ylabels("Source -> Target", fontsize=fs)
g.set_xlabels("Time", fontsize=fs)
g.set_axis_labels(size=fs, fontsize=fs)
for ax in g.axes.flat:
    ax.tick_params(axis="x", labelsize=fs)
    ax.tick_params(axis="y", labelsize=fs)

# Adjust layout for better readability
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Metastasis Times by Source-Target and Threshold", fontsize=fs)

# Save the grid of plots
plt.tight_layout()
plt.savefig(outfile)
plt.close()
