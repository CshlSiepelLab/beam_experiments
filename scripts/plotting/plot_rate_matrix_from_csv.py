#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = sys.argv[1]
outfile = sys.argv[2]


# file="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/model_selection_transition_matrices/model1.csv"
# outfile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/model_selection_transition_matrices/model1.pdf"


# Load the CSV file into a DataFrame
df = pd.read_csv(file, header=None)

state_labels = ["P", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"]

# Convert the DataFrame to a NumPy array
matrix = df.values

# Plot the matrix
num_rows, num_cols = matrix.shape
textsize = 18

# Set the diagonal values to NaN
np.fill_diagonal(matrix, np.nan)

# Plot the matrix with empty boxes for diagonal values
plt.imshow(matrix, cmap="viridis", interpolation="nearest", vmin=0, vmax=1)
for i in range(num_rows + 1):
    plt.hlines(i - 0.5, xmin=-0.5, xmax=num_cols - 0.5, color="white", linewidth=1)
for j in range(num_cols + 1):
    plt.vlines(j - 0.5, ymin=-0.5, ymax=num_rows - 0.5, color="white", linewidth=1)
plt.imshow(matrix, cmap="viridis", interpolation="nearest", vmin=0, vmax=1)
for i in range(num_rows + 1):
    plt.hlines(i - 0.5, xmin=-0.5, xmax=num_cols - 0.5, color="white", linewidth=1)
for j in range(num_cols + 1):
    plt.vlines(j - 0.5, ymin=-0.5, ymax=num_rows - 0.5, color="white", linewidth=1)

# Add index and column labels for the groups
plt.xticks(range(num_cols), state_labels, fontsize=textsize, rotation=90, ha="center")
plt.gca().xaxis.tick_top()
plt.yticks(range(num_rows), state_labels, fontsize=textsize)
plt.title("Recipient", fontsize=textsize)
plt.xlabel("", fontsize=textsize)
plt.ylabel("Source", fontsize=textsize)
plt.colorbar()
plt.tight_layout()
# plt.show()

plt.savefig(outfile)

plt.close()
