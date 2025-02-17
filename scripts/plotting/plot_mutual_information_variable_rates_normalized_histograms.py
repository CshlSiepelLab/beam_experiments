#!/usr/bin/env python3

import pandas as pd
import sys
import matplotlib.pyplot as plt


# input_file_path = sys.argv[1]
# output_file_path = sys.argv[2]

input_file_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_1_31_25_data_from_8_19_24/gtr_beam_mutual_information.csv"
output_file_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_1_31_25_data_from_8_19_24/gtr_beam_mutual_information_variable_rates.pdf"

# Load the data
data = pd.read_csv(input_file_path)

# Extract migration and mutation rates from simname
data[['mig_rate', 'mut_rate', 'run_id']] = data['simname'].str.extract(r'mig(\d+)_mut(\d+)_(\d+)')

# Convert migration rates to numeric
data['mig_rate'] = data['mig_rate'].apply(lambda x: float('1e-' + x))
data['mut_rate'] = data['mut_rate'].apply(lambda x: float('0.' + x))

# Get unique rates
mig_rates = sorted(data['mig_rate'].unique())
mut_rates = sorted(data['mut_rate'].unique())

# Create subplot grid
fig, axes = plt.subplots(len(mig_rates), len(mut_rates), figsize=(15, 15))
plt.subplots_adjust(hspace=0.4, wspace=0.3)

fs=22

# Plot histograms for each combination
for i, mig in enumerate(mig_rates):
    for j, mut in enumerate(mut_rates):
        subset = data[(data['mig_rate'] == mig) & (data['mut_rate'] == mut)]
        if not subset.empty:
            axes[i,j].hist(subset['mutual_information_normalized'], bins=20, range=(0,1), edgecolor='black', color='gray')
            # Only set the migration rate label for the leftmost column
            if j == 0:
                axes[i,j].set_ylabel(f'{mig}', fontsize=fs)
            # Only set the mutation rate label for the top row
            if i == 0:
                axes[i,j].set_title(f'{mut:.4f}', fontsize=fs)
            axes[i,j].tick_params(labelsize=fs/2)
            axes[i,j].tick_params(axis='x', rotation=0)
            axes[i,j].set_xticks([0, 0.25, 0.5, 0.75, 1.0])

# Add overall labels
fig.supylabel('Migration Rate', fontsize=fs, x=0.05)
fig.supxlabel('Mutual Information', fontsize=fs, y=0.05)
fig.suptitle('Mutation Rate', fontsize=fs, y=0.925)

plt.savefig(output_file_path)
plt.close()