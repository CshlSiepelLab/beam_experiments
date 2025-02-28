#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


file_path = sys.argv[1]
outfile = sys.argv[2]


# file_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/all_consensus_classifications.csv'
# outfile = "./test.pdf"


data = pd.read_csv(file_path, index_col=0)

# Calculate the percentage of True values in each column
percentages = data.mean() * 100

# Create a bar plot
plt.figure(figsize=(10, 6))
sns.barplot(x=percentages.index, y=percentages.values, palette="Greys")

plt.ylim(0,100)

plt.xlabel('Seeding topology', fontsize=22)
plt.ylabel('Percent of data', fontsize=22)
plt.xticks(rotation=0, fontsize=22)
plt.yticks(fontsize=22)

# Save the plot to a file
plt.tight_layout()
plt.savefig(outfile)
plt.close()

