#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path):
    return pd.read_csv(file_path, index_col=0)

def calculate_percentages(data, method_name):
    percentages = data.mean().reset_index()
    percentages.columns = ['Seeding topology', 'Value']
    percentages['Value'] *= 100
    percentages['method'] = method_name
    return percentages

def plot_percentages(percentages, outfile):
    # Rename the seeding topologies
    percentages['Seeding topology'] = percentages['Seeding topology'].replace({
        'met_to_met': 'Met to Met',
        'met_to_primary': 'Primary\nReseeding'
    })

    # Order the data to have the bars with MACH2 first and then BEAM
    percentages['method'] = pd.Categorical(percentages['method'], categories=['MACH2', 'BEAM'], ordered=True)

    # Create a bar plot
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Seeding topology', y='Value', hue='method', data=percentages, palette=["Blue", "Green"])

    plt.ylim(0, 100)
    plt.ylabel('Percent of data', fontsize=22)
    plt.xlabel('')
    plt.xticks(rotation=0, fontsize=22)
    plt.yticks(fontsize=22)
    plt.legend(title='Method', fontsize=18, title_fontsize=20, bbox_to_anchor=(1.05, 0.5), loc='upper left', borderaxespad=0., frameon=False)

    # Save the plot to a file
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def main():
    file_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/all_consensus_classifications.csv"
    file_path2 = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/mach2/5k/all_consensus_classifications.csv"
    outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/all_consensus_classifications.pdf"

    data = load_data(file_path)
    data2 = load_data(file_path2)

    percentages1 = calculate_percentages(data, 'BEAM')
    percentages2 = calculate_percentages(data2, 'MACH2')

    # Combine the percentages dataframes
    percentages = pd.concat([percentages1, percentages2])

    plot_percentages(percentages, outfile)

if __name__ == "__main__":
    main()
