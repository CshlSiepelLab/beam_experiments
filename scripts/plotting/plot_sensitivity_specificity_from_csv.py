#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def calculate_sensitivity_specificity(df, method):
    true_positive = np.sum((df['true'] == 'yes') & (df[method] == 'yes'))
    true_negative = np.sum((df['true'] == 'no') & (df[method] == 'no'))
    false_positive = np.sum((df['true'] == 'no') & (df[method] == 'yes'))
    false_negative = np.sum((df['true'] == 'yes') & (df[method] == 'no'))

    sensitivity = true_positive / (true_positive + false_negative)
    specificity = true_negative / (true_negative + false_positive)

    return sensitivity, specificity


# Load the CSV file
file_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/model_selection_2_or_3_parameters_reseeding_no_reseeding_12_7_24_variable_rates_data_8_19_24/all_classification_results.csv'

df = pd.read_csv(file_path)

# Classify beam Bayes factors as supporting reseeding or not
min_bf = 5 # 0 = support barely worth mentioning, 1.1 = positive support, 3 = strong support, 5 = overwhelming support
df['beam'] = df.apply(lambda row: 'yes' if row['beam_bf'] > min_bf and row['beam_bf'] > row['min_bf_diff'] else 'no', axis=1)

# Drop rows with NaN values
df = df.fillna('no')    # assumes those that did not converge are not supporting reseeding. This is a safe assumption since the null hypothesis is no reseeding.

methods = [name for name in df.columns if name not in ['sim_name', 'true', 'beam_bf', 'min_bf_diff', 'information']]
sensitivity = []
specificity = []

for method in methods:
    sens, spec = calculate_sensitivity_specificity(df, method)
    sensitivity.append(sens)
    specificity.append(spec)

# Plotting params
fs = 24
colors = ["#6fa8dc", "#e69138", "#c27ba0"]
methods = [method.capitalize() for method in methods]

# Plot two distributions of the beam_bf values on the same plot where each distribution is defined by the true classification
true_yes = df[df['true'] == 'yes']['beam_bf']
true_no = df[df['true'] == 'no']['beam_bf']

# Remove outliers more than 100
true_yes = true_yes[true_yes < 100]
true_no = true_no[true_no < 100]

plt.figure(figsize=(10, 6))
bins = 100
sns.kdeplot(true_yes, label='Reseeding\nground truth', color='blue', fill=True, alpha=0.5)
sns.kdeplot(true_no, label='No reseeding\nground truth', color='orange', fill=True, alpha=0.5)
plt.xlabel('Bayes Factor', fontsize=fs)
plt.ylabel('Frequency', fontsize=fs)
plt.legend(fontsize=fs, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tight_layout()
outfile_dist = file_path.replace('.csv', '_beam_bf_distribution.pdf')
plt.savefig(outfile_dist)
plt.close()

# Plotting sensitivity and specificity bar plots
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

axes[0].bar(methods, sensitivity, color=colors)
axes[0].set_ylim(0, 1)
axes[0].set_ylabel('Sensitivity', fontsize=fs)
axes[0].tick_params(axis='both', which='major', labelsize=fs)

axes[1].bar(methods, specificity, color=colors)
axes[1].set_ylim(0, 1)
axes[1].set_ylabel('Specificity', fontsize=fs)
axes[1].tick_params(axis='both', which='major', labelsize=fs)

plt.tight_layout()
outfile = file_path.replace('.csv', '_sensitivity_specificity.pdf')
plt.savefig(outfile)
plt.close()

# QC to report any sims needing more nested sampling runs
sims_needing_particles = df[abs(df['beam_bf']) < df['min_bf_diff']]['sim_name'].tolist()
print(f"{len(sims_needing_particles)} sims need more particles to decrease the minimum Bayes factor difference threshold.")
print(sims_needing_particles)
