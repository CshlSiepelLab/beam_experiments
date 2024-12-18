#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV file
# file_path = sys.argv[1]

file_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/beam_information_content.csv'

outfile = file_path.replace('.csv', '.pdf')

df = pd.read_csv(file_path)

df['mig'] = df['sim_name'].str.extract(r'mig(\d+)')[0]
df['mut'] = df['sim_name'].str.extract(r'mut(\d+)')[0]

df['mut'] = [f"0.{val}" for val in df['mut']]
df['mut'] = df['mut'].astype(float)

df['mig'] = [f"1e-{val}" for val in df['mig']]


unique_migs = df['mig'].unique()
unique_muts = df['mut'].unique()

df = df.sort_values(by=['mut'], ascending=True)
df = df.sort_values(by=['mig'], ascending=False)

# Plotting
fs = 26
plt.figure(figsize=(12, 8))

sns.boxplot(x='mut', y='information_content', hue='mig', data=df, showcaps=True, showfliers=False)
sns.stripplot(x='mut', y='information_content', hue='mig', color="black", data=df, dodge=True, alpha=0.5, legend=False)

plt.xlabel('Mutation rate', fontsize=fs)
plt.ylabel('Information', fontsize=fs)
plt.legend(title='Migration rate', title_fontsize=fs, fontsize=fs, bbox_to_anchor=(1.05, 0.7), loc='upper left', frameon=False)
plt.tick_params(axis='both', which='major', labelsize=fs)

plt.tight_layout()
plt.savefig(outfile)
plt.close()
