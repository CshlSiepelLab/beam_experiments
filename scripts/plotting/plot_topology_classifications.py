#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

infile = sys.argv[1]

# testing
# infile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/compare_migration_counts/topology_classifications_all.csv"

# Load the data into a DataFrame
data = pd.read_csv(infile)

# Calculate the percentage of True values for each category and group
summary = {
    'group': ['machina', 'machina', 'metient', 'metient', 'metastabayes', 'metastabayes'],
    'category': ['met_to_met', 'reseeding', 'met_to_met', 'reseeding', 'met_to_met', 'reseeding'],
    'percentage': [
        data['machina_met_to_met'].mean() * 100,
        data['machina_reseeding'].mean() * 100,
        data['metient_met_to_met'].mean() * 100,
        data['metient_reseeding'].mean() * 100,
        data['metastabayes_met_to_met'].mean() * 100,
        data['metastabayes_reseeding'].mean() * 100
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
    'metastabayes': 'Beast'
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
