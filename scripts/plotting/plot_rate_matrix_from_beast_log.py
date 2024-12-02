#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import re
import matplotlib.pyplot as plt

# input_file = sys.argv[1]
# primary_tissue = sys.argv[2]

input_file = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1467/CP01/combined.log'
primary_tissue = 'PRL'

outfile = input_file.replace(".log", "_tissue_rate_matrix.pdf")

df = pd.read_csv(input_file, sep='\t', comment='#')

tissue_rate_col_names = [name for name in df.columns if name.startswith('geoSubstModelLogger')]
tissues = list(set([tissue for name in tissue_rate_col_names for tissue in name.replace("geoSubstModelLogger.relGeoRate_", "").split("_") ]))
tissues = [primary_tissue] + sorted([tis for tis in tissues if tis != primary_tissue])

# Calculate average rates
rates = {}
# total_sum = 0  

for source in tissues:
    sum_rates = 0
    for recipient in tissues:
        if source == recipient:
            continue
        rates[f"{source}_{recipient}"] = df[f'geoSubstModelLogger.relGeoRate_{source}_{recipient}'].mean()

#         sum_rates += rates[f"{source}_{recipient}"]
    
#     rates[f"{source}_{source}"] = -sum_rates
#     total_sum += sum_rates

# # Normalize rates by dividing by total_sum
# for key in rates:
#     rates[key] /= total_sum

# Create a matrix for the rates
num_tissues = len(tissues)
rate_matrix = np.zeros((num_tissues, num_tissues))

for (source_recipient), rate in rates.items():
    source, recipient = source_recipient.split('_')
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

# plot the distribution histogram for geoSubstModelLogger.relGeoRate columns for those rates with the recipient as primary_tissue
primary_tissue_rates = [df[f'geoSubstModelLogger.relGeoRate_{source}_{primary_tissue}'] for source in tissues if source != primary_tissue]

for source in tissues:
    if source == primary_tissue:
        continue
    primary_tissue_rate = df[f'geoSubstModelLogger.relGeoRate_{source}_{primary_tissue}']

    plt.figure(figsize=(8, 6))
    sns.histplot(primary_tissue_rate, bins=30, kde=False, color='grey')
    plt.xlabel('Rate', fontsize=fs)
    plt.ylabel('Frequency', fontsize=fs)
    plt.title(f'{source} to {primary_tissue}', fontsize=fs)
    plt.xticks(fontsize=ts)
    plt.yticks(fontsize=ts)
    plt.tight_layout()

    histogram_outfile = input_file.replace(".log", f"_{source}_to_{primary_tissue}_rate_distribution.pdf")
    plt.savefig(histogram_outfile)

    plt.close()