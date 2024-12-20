#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

infile = sys.argv[1]

# # testing
# infile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/compare_migration_counts/migration_counts_all.csv"

# Load the data into a DataFrame
data = pd.read_csv(infile)

# # Discard any rows with metastabayes_migrations as 0
# data = data[data['metastabayes_migrations'] != 0]

# Melt the DataFrame to long format
melted_data = pd.melt(data, id_vars=['name'], 
                      value_vars=['parsimony_migrations', 'parsimony_comigrations',
                                'machina_migrations', 'machina_comigrations', 
                                'metient_migrations', 'metient_comigrations', 
                                'beam_migrations', 'beam_comigrations'],
                      var_name='category', value_name='count')

# Extract the group and type from the 'category' column
melted_data['group'] = melted_data['category'].apply(lambda x: x.split('_')[0])
melted_data['type'] = melted_data['category'].apply(lambda x: x.split('_')[1])

# Rename methods
melted_data['group'] = melted_data['group'].replace({'parsimony': 'Parsimony', 'machina': 'Machina', 'metient': 'Metient', 'beam': 'Beam'})
melted_data['type'] = melted_data['type'].replace({'migrations': 'Migrations', 'comigrations': 'Co-migrations'})

# Remove rows where 'group' is 'Parsimony'
melted_data = melted_data[melted_data['group'] != 'Parsimony']

DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

# Create the boxplot
plt.figure(figsize=(12, 8))

sns.boxplot(x='group', y='count', hue='type', data=melted_data, palette=DEFAULT_COLORS[2:], showfliers=False, linewidth=2)
# sns.stripplot(x='group', y='count', hue='type', data=melted_data, dodge=True, color='grey', alpha=0.5, legend=False)

# Set the labels and title
plt.xlabel('', fontsize=28)
plt.ylabel('Count', fontsize=28)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)

# Show the plot
plt.legend(title='', bbox_to_anchor=(1.0, 0.6), loc='upper left', fontsize=24, frameon=False)
plt.tight_layout()

plt.savefig(infile.replace('.csv', '.pdf'))
