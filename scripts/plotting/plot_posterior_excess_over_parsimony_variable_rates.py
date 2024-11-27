#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
csv_file = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/metastabayes/all_posterior_expected_migration_counts_over_parsimony.csv'

data = pd.read_csv(csv_file)

# Calculate the difference between the two columns
data['difference'] = data['beast_migration_count'] - data['parsimony_migration_count']
data['original_name'] = data['name']
data[['mig', 'mut', 'name']] = data['name'].str.split('_', expand=True)

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
grouped = data.groupby('original_name')['difference'].mean().reset_index()
plt.figure(figsize=(10, 6))
plt.hist(grouped['difference'], bins=100, edgecolor='black', color='grey')
plt.xlabel('Posterior expected\n[BEAM migration count - Parsimony migration count]', fontsize=fontsize)
plt.ylabel('Number of datasets', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.xlim(0, 25)
# plt.ylim(0, 35)
# plt.xticks(np.arange(0, 7.5, 1), fontsize=fontsize)
# plt.yticks(np.arange(0, 40, 5), fontsize=fontsize)
plt.tight_layout()
outfile = csv_file.replace('.csv', '_difference_means_histogram.pdf')
plt.savefig(outfile)
plt.close()

# Make one plot for each pair of mig and mut from the already grouped data
grouped[['mig', 'mut', 'name']] = grouped['original_name'].str.split('_', expand=True)
for mig in grouped['mig'].unique():
    for mut in grouped['mut'].unique():
        subset = grouped[(grouped['mig'] == mig) & (grouped['mut'] == mut)]
        plt.figure(figsize=(6, 6))
        plt.hist(subset['difference'], bins=15, edgecolor='black', color='grey')
        plt.xlabel('Posterior expected\n[BEAM migration count - Parsimony migration count]', fontsize=fontsize)
        plt.ylabel('Number of datasets', fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlim(0, 22)
        # plt.ylim(0, subset['difference'].max() + 1)
        plt.xticks(np.arange(0, 23, 5), fontsize=fontsize)
        max_y = plt.ylim()[1]
        if max_y < 10:
            inc = 1
        else:
            inc = 5
        plt.yticks(np.arange(0, max_y + 1, inc), fontsize=fontsize)
        plt.tight_layout()
        outfile = csv_file.replace('.csv', f'_difference_means_histogram_mig_{mig}_mut_{mut}.pdf')
        plt.savefig(outfile)
        plt.close()