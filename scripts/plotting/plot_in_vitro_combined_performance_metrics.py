#!/usr/bin/env python3

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# infile = sys.argv[1]
# outprefix = sys.argv[2]

infile = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25/combined_performance_metrics.csv"
prefix = "subset_38_phylogeneticallyInformative_"
outprefix = os.path.join(os.path.dirname(infile), prefix + os.path.basename(infile).replace(".csv", ""))

df = pd.read_csv(infile)

# Only keep CPs solved by all methods for each scenario
num_methods = df['method'].nunique()
df = df.groupby(['scenario', 'cp']).filter(lambda x: x['method'].nunique() == num_methods)

# Set scenario col as int
df['scenario'] = df['scenario'].astype(int)

# Get the unique scenarios to plot seperately
scenarios = df['scenario'].unique()

# Set CPs to keep based on each scenario
# cps = {
#     1: ['CP020', 'CP021', 'CP032', 'CP035', 'CP027', 'CP013', 'CP054', 'CP018', 'CP038'],
#     2: ['CP021', 'CP035', 'CP054', 'CP017', 'CP023'],
#     3: ['CP021', 'CP035', 'CP054']
# }
cps = ['CP011', 'CP014', 'CP019', 'CP016', 'CP015', 'CP013', 'CP012', 'CP024', 'CP018', 'CP020', 'CP017', 'CP021', 'CP023', 'CP022', 'CP025', 'CP026', 'CP036', 'CP030', 'CP029', 'CP035', 'CP093', 'CP045', 'CP038', 'CP034', 'CP027', 'CP031', 'CP028', 'CP040', 'CP033', 'CP065', 'CP061', 'CP054', 'CP037', 'CP032', 'CP190', 'CP121', 'CP055']

fs = 22

for scenario in scenarios:

    # Specify cps to keep
    # df = df[df['cp'].isin(cps[scenario])]
    df = df[df['cp'].isin(cps)]

    df_scenario = df[df['scenario'] == scenario]

    # Always have mach2 before beam
    order = ['mach2', 'beam']
    
    # Downsample data thresholds
    thresholds = list(np.arange(0, 1.1, 0.1)) + [0.95, 0.99]
    df_scenario = df_scenario[df_scenario['threshold'].isin(thresholds)]
    df_scenario = df_scenario.sort_values(by='method')
    
    # Accuracy plot
    plt.figure(figsize=(15, 8))
    sns.boxplot(data=df_scenario, x='threshold', y='accuracy', hue='method', hue_order=order, showfliers=False, boxprops=dict(edgecolor="black", alpha=0.9),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        medianprops=dict(color="black"))
    sns.stripplot(data=df_scenario, x='threshold', y='accuracy', hue='method', hue_order=order, palette=["black", "black"], size=5, dodge=True, alpha=0.5, legend=False)
    plt.title('', fontsize=fs)
    plt.xlabel('Consensus graph threshold', fontsize=fs)
    plt.ylabel('Accuracy', fontsize=fs)
    plt.xticks(rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fs, title_fontsize=fs)
    plt.tight_layout()
    plt.savefig(f"{outprefix}_accuracy_scenario{scenario}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # F1 plot
    plt.figure(figsize=(15, 8))
    sns.boxplot(data=df_scenario, x='threshold', y='f1', hue='method', hue_order=order, showfliers=False, boxprops=dict(edgecolor="black", alpha=0.9),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        medianprops=dict(color="black"))
    sns.stripplot(data=df_scenario, x='threshold', y='f1', hue="method", hue_order=order, palette=["black", "black"], size=5, alpha=0.5, dodge=True, legend=False)
    plt.title('', fontsize=fs)
    plt.xlabel('Consensus graph threshold', fontsize=fs)
    plt.ylabel('F1', fontsize=fs)
    plt.xticks(rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fs, title_fontsize=fs)
    plt.tight_layout()
    plt.savefig(f"{outprefix}_f1_scenario{scenario}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    df_scenario = df[df['scenario'] == scenario]
    df_scenario = df_scenario.sort_values(by='method')
    
    # Plot precision recall curves
    plt.figure(figsize=(10, 10))
    
    # # Calculate mean precision and recall for each method
    # for method in df_scenario['method'].unique():
    #     method_data = df_scenario[df_scenario['method'] == method]
    #     # Calculate means
    #     mean_precision = method_data.groupby('threshold')['precision'].mean()
    #     mean_recall = method_data.groupby('threshold')['accuracy'].mean()
        
    #     # Create sorted arrays for plotting
    #     sorted_data = pd.DataFrame({
    #         'recall': mean_recall,
    #         'precision': mean_precision
    #     }).sort_values('recall')
        
    #     plt.plot(sorted_data['recall'], sorted_data['precision'], 
    #             marker='o', label=method, linewidth=2, markersize=0)
    
    # plt.xlabel('Recall', fontsize=fs)
    # plt.ylabel('Precision', fontsize=fs)
    # plt.xticks(fontsize=fs)
    # plt.yticks(fontsize=fs)
    # plt.xlim(0,1)
    # plt.ylim(0,1)
    # plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fs, title_fontsize=fs)
    # plt.grid(True, linestyle='--', alpha=0.7)
    # plt.tight_layout()
    # plt.savefig(f"{outprefix}_pr_curve_scenario{scenario}.pdf", dpi=300, bbox_inches='tight')
    # plt.close()
    
