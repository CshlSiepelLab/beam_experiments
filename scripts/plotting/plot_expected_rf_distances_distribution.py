#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt


input_csv_path = sys.argv[1]
output_plot_path = sys.argv[2]

# Read the CSV file
data = pd.read_csv(input_csv_path)

# Extract the 'expected_rf_distance' column
distances = data["expected_rf_distance"]

# Plot the histogram
plt.figure(figsize=(8, 6))

fs = 22

plt.hist(distances, bins=40, color="grey", edgecolor="black", alpha=1.0)
plt.xlabel(
    "E[Normalized RF distance between BEAM\nposterior sample and LAML tree]",
    fontsize=fs,
)
plt.ylabel("Number of CPs", fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlim(0, 1)

# Improve layout for publication
plt.tight_layout()
plt.savefig(output_plot_path)
plt.close()
