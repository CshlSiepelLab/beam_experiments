#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import re
import matplotlib.pyplot as plt

# input_file = sys.argv[1]
# primary_tissue = sys.argv[2]

input_file = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_1_20_25_uniform_50cells_50sites_data_7_24_24/beam_gtr/mS_12234/combined.log'
primary_tissue = 'P'

outfile = input_file.replace(".log", "_tissue_rate_matrix.pdf")

df = pd.read_csv(input_file, sep='\t', comment='#')

tissue_rate_col_names = [name for name in df.columns if name.startswith('tissueSubstModelLogger')]
tissues = list(set([tissue for name in tissue_rate_col_names for tissue in name.replace("tissueSubstModelLogger.relGeoRate_", "").split("_") ]))
tissues = [primary_tissue] + sorted([tis for tis in tissues if tis != primary_tissue])

# Get the rates (note they are already logged with respect to the equillibrium frequencies of each tissue for normalization to an expectation of one substitution per unit time)
rates = {}
sample = np.random.randint(0, 1000) # if using a single sample, then draw a random sample from the end of the posterior samples
print(sample)
# Create a matrix for the rates
num_tissues = len(tissues)
rate_matrix = np.zeros((num_tissues, num_tissues))

for source in tissues:
    for recipient in tissues:
        if source == recipient:
            continue
        rate = df[f'tissueSubstModelLogger.relGeoRate_{source}_{recipient}'].iloc[-sample] # to use a single sample
        # rate = df[f'tissueSubstModelLogger.relGeoRate_{source}_{recipient}'].mean() # to use the mean value
        i = tissues.index(source)
        j = tissues.index(recipient)
        rate_matrix[i, j] = rate

# Plot the rate matrix
fs = 28
ts = 18
plt.figure(figsize=(10, 8))
heatmap = sns.heatmap(rate_matrix, xticklabels=tissues, yticklabels=tissues, annot=True, cmap='YlOrRd', annot_kws={"size": ts}, cbar_kws={'label': 'Rate'})
heatmap.figure.axes[-1].yaxis.label.set_size(fs)
heatmap.figure.axes[-1].tick_params(labelsize=ts)
plt.xlabel('Recipient', fontsize=fs)
plt.ylabel('Source', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.tight_layout()

plt.savefig(outfile)

plt.close()

# # plot the distribution histogram for tissueSubstModelLogger.relGeoRate columns for those rates with the recipient as primary_tissue
# primary_tissue_rates = [df[f'tissueSubstModelLogger.relGeoRate_{source}_{primary_tissue}'] for source in tissues if source != primary_tissue]

# for source in tissues:
#     if source == primary_tissue:
#         continue
#     primary_tissue_rate = df[f'tissueSubstModelLogger.relGeoRate_{source}_{primary_tissue}']

#     plt.figure(figsize=(8, 6))
#     sns.histplot(primary_tissue_rate, bins=30, kde=False, color='grey')
#     plt.xlabel('Rate', fontsize=fs)
#     plt.ylabel('Frequency', fontsize=fs)
#     plt.title(f'{source} to {primary_tissue}', fontsize=fs)
#     plt.xticks(fontsize=ts)
#     plt.yticks(fontsize=ts)
#     plt.tight_layout()

#     histogram_outfile = input_file.replace(".log", f"_{source}_to_{primary_tissue}_rate_distribution.pdf")
#     plt.savefig(histogram_outfile)

#     plt.close()