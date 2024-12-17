#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

# Read in the CSV files
performance_threshold_csv = sys.argv[1]
information_csv = sys.argv[2]
outfile = sys.argv[3]

# performance_threshold_csv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/bin_information_performance/concat_all_threshold_stats.csv"
# information_csv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/beam_information_content.csv"
# outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/bin_information_performance/bin_information_precision_recall.pdf"

df1 = pd.read_csv(performance_threshold_csv)
df2 = pd.read_csv(information_csv)

# Merge the dataframes on sim_name
df = pd.merge(df1, df2, left_on='sim', right_on='sim_name')

# Bin by information content
# Determine the bin size to have roughly 50 sims per bin
desired_bin_size = 80
total_sims = df['sim_name'].nunique()
num_bins = total_sims // desired_bin_size

# Create bins based on the desired number of bins
df['information_bin'] = pd.qcut(df['information_content'], num_bins, duplicates='drop')

# Obtain a dict mapping of information_bin value to the number of unique sim_name values for that bin
bin_count = df.groupby('information_bin')['sim_name'].nunique().to_dict()

# Group by 'information_bin' and 'threshold', then calculate the average precision and recall for each bin and threshold
grouped = df.groupby(['information_bin', 'Threshold'])[['precision', 'recall']].mean().reset_index()

grouped = grouped.sort_values(by=['information_bin', 'Threshold'])

# Drop 1.0 threshold
grouped = grouped[grouped['Threshold'] != 1.0]

# Plot the precision and recall curves for each threshold
fs = 32
plt.figure(figsize=(12, 8))
for bin in grouped['information_bin'].unique():
    subset = grouped[grouped['information_bin'] == bin]
    plt.plot(subset['recall'], subset['precision'], label=f"{bin.left} - {bin.right}\n(N = {bin_count[bin]})", linewidth=2.5)

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel('Recall', fontsize=fs)
plt.ylabel('Precision', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ts=22
plt.legend(title='Information bin', bbox_to_anchor=(1.05, 0.8), loc='upper left', fontsize=ts, edgecolor='none', title_fontsize=ts)
plt.tight_layout()
plt.savefig(outfile)
plt.close()