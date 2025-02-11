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
    fs=18

    colors = {
        'BEAM': 'lightgrey',
        'Random': 'darkgrey',
        'Parsimony': 'purple',
        'MACHINA': 'red',
        'Metient': 'green',
        'Consensus': 'blue',
        'MACH2' : 'navy',
        'FitchCount' : 'orange'
    }
    colors = [colors[col] for col in data_f1.columns]

    sns.boxplot(data=data_f1, orient='v', showfliers=False, palette=colors, boxprops=dict(edgecolor="black", alpha=0.9), whiskerprops=dict(color="black"), capprops=dict(color="black"), medianprops=dict(color="black"))
    sns.stripplot(data=data_f1, orient='v', color='black', alpha=0.5)
    plt.ylabel('F1 score', fontsize=fs)
    plt.xticks(rotation=45, fontsize=fs)
    plt.yticks(fontsize=fs)
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