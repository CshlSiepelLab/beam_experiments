#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

infile = sys.argv[1]

# # testing
# infile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/compare_migration_counts/topology_classifications_all.csv"

# Load the data into a DataFrame
data = pd.read_csv(infile)

# Calculate the percentage of True values for each category and group
summary = {
    'group': ['machina', 'machina', 'metient', 'metient', 'beam', 'beam'],
    'category': ['met_to_met', 'reseeding', 'met_to_met', 'reseeding', 'met_to_met', 'reseeding'],
    'percentage': [
        data['machina_met_to_met'].mean(skipna=True) * 100,
        data['machina_reseeding'].mean(skipna=True) * 100,
        data['metient_met_to_met'].mean(skipna=True) * 100,
        data['metient_reseeding'].mean(skipna=True) * 100,
        data['beam_met_to_met'].mean(skipna=True) * 100,
        data['beam_reseeding'].mean(skipna=True) * 100
    ]
}

summary_df = pd.DataFrame(summary)

# Rename the categories
summary_df['category'] = summary_df['category'].replace({
    'met_to_met': 'Met to Met',
    'reseeding': 'Reseeding'
})

summary_df['group'] = summary_df['group'].replace({
    'machina': 'Machina',
    'metient': 'Metient',
    'beam': 'Beam'
})

DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

# Create the bar plot
plt.figure(figsize=(10, 6))
sns.barplot(x='percentage', y='group', hue='category', data=summary_df, width = 0.75, palette=DEFAULT_COLORS[2:], orient='h')

plt.xlim(0, 100)

# Set the labels and title
plt.xlabel('Percent of topology across all CPs', fontsize=24)
plt.ylabel('', fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

# Show the plot
plt.legend(title='', bbox_to_anchor=(1.0, 0.6), loc='upper left', frameon=False, fontsize=20)
plt.tight_layout()

plt.savefig(infile.replace('.csv', '.pdf'))

plt.close()


# # Create Venn diagrams for each method to see if clonal populations have either met to met or reseeding OR they tend to have both or none
# methods = ['machina', 'metient', 'beam']
# method_names = ['Machina', 'Metient', 'Beam']

# fig, axes = plt.subplots(3, 1, figsize=(10, 10))

# for i, (method, method_name) in enumerate(zip(methods, method_names)):
#     met_to_met = data[f'{method}_met_to_met']
#     reseeding = data[f'{method}_reseeding']
    
#     # Calculate the sets for Venn diagram
#     set_met_to_met = set(data.index[met_to_met])
#     set_reseeding = set(data.index[reseeding])
    
#     plt.sca(axes[i])
#     venn2([set_met_to_met, set_reseeding], ('Met to Met', 'Reseeding'), set_colors=(DEFAULT_COLORS[2], DEFAULT_COLORS[3]))
#     for patch in axes[i].patches:
#         if patch is not None:
#             patch.set_alpha(0.5)
#             patch.set_edgecolor('black')
#             patch.set_linewidth(1.5)
#     for text in axes[i].texts:
#         text.set_fontsize(18)
#     axes[i].set_aspect(0.5)
#     axes[i].set_title(f'{method_name}', fontsize=24)

# plt.tight_layout()
# plt.savefig(infile.replace('.csv', '_venn_diagrams.pdf'))
# plt.close()
