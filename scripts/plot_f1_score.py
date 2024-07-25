#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

def plot_f1_score(csv_path, output_dir):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_path)
    
    # Filter columns with "f1" in the name
    f1_columns = [col for col in df.columns if 'f1' in col]

    data_f1 = df[f1_columns]
    # Remove "_f1" from all column names
    data_f1.columns = [col.replace('_f1', '') for col in data_f1.columns]

    # Create a boxplot for all f1 columns
    plt.figure()
    fs=16
    sns.boxplot(data=data_f1, orient='v')
    sns.stripplot(data=data_f1, orient='v', color='black', alpha=0.5)
    plt.xlabel('Method', fontsize=fs)
    plt.ylabel('F1 score \nof migration graph', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.xticks(rotation=45, fontsize=fs)
    plt.ylim(-0.05,1.05)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/f1_score.pdf')
    plt.close()


# user inputs
csv_path = sys.argv[1]
output_dir = sys.argv[2]

# # Example usage
# csv_path = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/precision_recall_curve/metrics.csv'
# output_dir = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/precision_recall_curve'

# plot f1
plot_f1_score(csv_path, output_dir)