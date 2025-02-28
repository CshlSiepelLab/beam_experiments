#!/usr/bin/env python3

import pandas as pd
from scipy.stats import ttest_rel, ttest_ind

# Load the CSV file
file_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_2_20_25_uniform_50cells_50sites_data_7_24_24/precision_recall_curve/metrics.csv'
data = pd.read_csv(file_path)

# Drop rows with NaN values in the relevant columns
cols = ['Metient_f1', 'BEAM_f1']
data = data.dropna(subset=cols)

# Subset the data
compare_f1 = data[cols[0]]
beam_f1 = data[cols[1]]

# Perform the paired t-test
t_stat, p_value = ttest_rel(compare_f1, beam_f1)

# Print the results
print(f'T-statistic: {t_stat}')
print(f'P-value: {p_value}')