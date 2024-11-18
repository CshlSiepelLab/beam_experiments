#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
csv_file = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/metastabayes/all_posterior_expected_migration_counts_over_parsimony.csv'

data = pd.read_csv(csv_file)

# Calculate the difference between the two columns
data['difference'] = data['beast_migration_count'] - data['parsimony_migration_count']

# Plot the distribution of the differences
fontsize = 24

plt.figure(figsize=(10, 6))
plt.hist(data['difference'], bins=100, edgecolor='black', color = 'grey')
plt.xlabel('All raw (BEAM - parsimony)\nmigration count differences', fontsize=fontsize)
plt.ylabel('Number of posterior samples', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.tight_layout()
outfile = csv_file.replace('.csv', '_all_differences_histogram.pdf')
plt.savefig(outfile)
plt.close()

# Plot the average difference grouped by sim_name, so compute the mean difference for each group and plot those values
grouped = data.groupby('name').mean()
plt.figure(figsize=(10, 6))
plt.hist(grouped['difference'], bins=100, edgecolor='black', color='grey')
plt.xlabel('Posterior expected\n[BEAM migration count - Parsimony migration count]', fontsize=fontsize)
plt.ylabel('Number of datasets', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.xlim(0, 7)
plt.ylim(0, 35)
plt.xticks(np.arange(0, 7.5, 1), fontsize=fontsize)
plt.yticks(np.arange(0, 40, 5), fontsize=fontsize)
plt.tight_layout()
outfile = csv_file.replace('.csv', '_difference_means_histogram.pdf')
plt.savefig(outfile)
plt.close()
