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
    data["num_mach2_solutions"],
    bins=range(1, 12),
    edgecolor="black",
    color="gray",
    align="left",
)
plt.xticks(range(1, 12))
fs = 22
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel("Number of MACH2 solutions", fontsize=fs)
plt.ylabel("Count", fontsize=fs)
plt.tight_layout()

# Save the plot to a PDF file
plt.savefig(output_file_path)
plt.close()
