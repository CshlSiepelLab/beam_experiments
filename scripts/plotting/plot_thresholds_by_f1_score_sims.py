#!/usr/bin/env python3

import sys
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

threshold_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_12_31_24_uniform_50cells_50sites_data_7_24_24/precision_recall_curve/all_threshold_stats.csv"
# threshold_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/all_threshold_stats.csv"

outdir = os.path.dirname(threshold_file)

df = pd.read_csv(threshold_file)

df['f1'] = 2 * (df['precision'] * df['recall']) / (df['precision'] + df['recall'])

df['Threshold'] = df['Threshold'].round(2)

# Group by threshold and calculate the average F1 score per threshold
average_f1_per_threshold = df.groupby('Threshold')['f1'].mean().reset_index()

# Rank the thresholds by their average F1 score
ranked_thresholds = average_f1_per_threshold.sort_values(by='f1', ascending=False)

# Find the threshold with the highest average F1 score
highest_f1_threshold = ranked_thresholds.iloc[0]

print("Threshold with the highest average F1 score:")
print(highest_f1_threshold)

plt.figure(figsize=(12, 6))
sns.scatterplot(x="Threshold", y="f1", data=df, size = 0.5, alpha = 0.25, color="grey", legend=False)
sns.lineplot(x="Threshold", y="f1", data=df, errorbar=None, estimator=np.median)
plt.xlabel('Posterior threshold', fontsize=24)
plt.ylabel('F1 Score', fontsize=24)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig(f"{outdir}/thresholds_f1_lineplot.pdf")