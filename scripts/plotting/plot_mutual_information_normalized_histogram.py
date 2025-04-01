#!/usr/bin/env python3

import pandas as pd
import sys
import matplotlib.pyplot as plt


input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# Load the data
data = pd.read_csv(input_file_path)

# Plot the histogram
plt.hist(
    data["mutual_information_normalized"],
    bins=100,
    range=(0, 1),
    edgecolor="black",
    color="gray",
)
plt.xlabel("Mutual information normalized")
plt.ylabel("Count")
fs = 22
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel("Mutual information normalized", fontsize=fs)
plt.ylabel("Count", fontsize=fs)
plt.tight_layout()

# Save the plot to a PDF file
plt.savefig(output_file_path)
plt.close()
