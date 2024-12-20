#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

def plot_f1_score(csv_path, output_dir):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_path)

    df['seeding_classification'] = [name.split("_")[0] for name in df['sim']]
    
    # Filter columns with "f1" in the name
    f1_columns = [col for col in df.columns if 'f1' in col and 'PathFinder' not in col]

    data_f1 = df[f1_columns]
    # Remove "_f1" from all column names
    data_f1.columns = [col.replace('_f1', '') for col in data_f1.columns]

    # Add seeding classification to the data_f1 DataFrame
    data_f1['seeding_classification'] = df['seeding_classification']

    # Create a boxplot for each seeding classification
    # unique_classes = data_f1['seeding_classification'].unique()
    unique_classes = ['mS', 'pS', 'pM', 'pR']
    num_classes = len(unique_classes)
    
    fig, axes = plt.subplots(1, num_classes, figsize=(5 * num_classes, 6), sharey=False)
    fs = 22

    for i, seeding_class in enumerate(unique_classes):
        class_data = data_f1[data_f1['seeding_classification'] == seeding_class].drop(columns=['seeding_classification'])
        sns.boxplot(data=class_data, orient='v', ax=axes[i])
        sns.stripplot(data=class_data, orient='v', color='black', alpha=0.5, ax=axes[i])
        axes[i].set_title(seeding_class, fontsize=fs)
        axes[i].set_ylabel('F1 score', fontsize=fs)
        axes[i].tick_params(axis='x', labelrotation=45, labelsize=fs)
        axes[i].tick_params(axis='y', labelsize=fs)
        axes[i].set_ylim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/f1_score_by_seeding_topology.pdf')
    plt.close()


# user inputs
csv_path = sys.argv[1]
output_dir = sys.argv[2]

# # Example usage
# csv_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/precision_recall_curve/metrics.csv'
# output_dir = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/precision_recall_curve'

# plot f1
plot_f1_score(csv_path, output_dir)