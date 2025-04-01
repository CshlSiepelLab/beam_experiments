#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# infile = sys.argv[1]

# testing
infile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/consensus_edge_probabilities_by_topology/consensus_graph_probabilities_by_edge_topology.csv"

outfile = infile.replace(".csv", ".pdf")

# Load the data into a DataFrame
data = pd.read_csv(infile)

# Melt the DataFrame to long format for seaborn
melted_data = data.melt(
    id_vars=["mmus", "cp"],
    value_vars=[
        "primary_seeding_avg_prob",
        "met_to_met_avg_prob",
        "primary_reseeding_avg_prob",
    ],
    var_name="Probability_Type",
    value_name="Probability",
)

# Rename the values in the 'Probability_Type' column
melted_data["Probability_Type"] = melted_data["Probability_Type"].replace(
    {
        "primary_seeding_avg_prob": "Primary to met",
        "met_to_met_avg_prob": "Met to met",
        "primary_reseeding_avg_prob": "Primary reseeding",
    }
)

DEFAULT_COLORS = ["grey", "#6fa8dc", "#e69138"]
fs = 24

# Create the boxplot
plt.figure(figsize=(8, 8))
sns.boxplot(
    x="Probability_Type",
    y="Probability",
    data=melted_data,
    showcaps=True,
    boxprops=dict(facecolor="none", edgecolor="black"),
    whiskerprops=dict(color="black"),
    medianprops=dict(color="black"),
)
sns.stripplot(
    x="Probability_Type",
    y="Probability",
    data=melted_data,
    palette=DEFAULT_COLORS,
    alpha=0.5,
    jitter=True,
    dodge=False,
    size=8,
)

# Set plot labels and title
plt.xlabel("")
plt.ylabel("Mean probability per edge\ntype in consensus graph", fontsize=fs)

# Set tick labels font size
plt.xticks(rotation=45, fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)

plt.tight_layout()
plt.savefig(outfile, bbox_inches="tight")
plt.close()
