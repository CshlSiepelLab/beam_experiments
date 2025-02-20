#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt

# File paths
file1 = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3/gtr_beam_mutual_information.csv'
file2 = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3/marginal_likelihoods.csv'
outfile = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3/correlation_mutual_info_bayes_factor.pdf'

# Read data
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Prepare and merge dataframes
df1['name'] = df1['mouse_cp']
df2['name'] = df2['name'].str.replace('/', '_')

df1 = df1[['name','mutual_information_normalized']]
df2 = df2[['name', 'bf(gtr-random)']]

# Merge dataframes on 'sim_name'
merged_df = pd.merge(df1, df2, on='name')

# Extract relevant columns
mutual_info = merged_df['mutual_information_normalized']
bayes_factor = merged_df['bf(gtr-random)']

# Plotting
plt.figure(figsize=(10, 6))
sns.scatterplot(x=mutual_info, y=bayes_factor, color='grey')
# sns.regplot(x=mutual_info, y=bayes_factor, scatter=False, color='red')

# Calculate Spearman correlation
spearman_corr, _ = spearmanr(mutual_info, bayes_factor)

# # Linear regression
# X = mutual_info.values.reshape(-1, 1)
# y = bayes_factor.values
# reg = LinearRegression().fit(X, y)
# r_squared = reg.score(X, y)

plt.xlabel('Mutual Information')
plt.ylabel('Bayes Factor')
plt.axhline(y=5, color='blue', linestyle='--')
# plt.title(f'y = {reg.coef_[0]:.2f}x + {reg.intercept_:.2f}\nR^2: {r_squared:.2f}\nSpearman Correlation: {spearman_corr:.2f}')
plt.tight_layout()

plt.savefig(outfile)
plt.close()
