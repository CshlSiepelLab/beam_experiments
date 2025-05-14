#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# infile = sys.argv[1]
# outprefix = sys.argv[2]

infile = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25/combined_performance_metrics.csv"
outprefix = infile.replace(".csv", "")

df = pd.read_csv(infile)

# Only keep CPs solved by all methods for each scenario
num_methods = df['method'].nunique()
df = df.groupby(['scenario', 'cp']).filter(lambda x: x['method'].nunique() == num_methods)

# Only keep CPs with phylogenetic information estimated seperately
cps = ['CP001', 'CP002', 'CP003', 'CP004', 'CP005', 'CP006', 'CP007', 'CP008', 'CP009', 'CP010', 'CP011', 'CP012', 'CP013', 'CP014', 'CP015', 'CP016', 'CP017', 'CP018', 'CP019', 'CP020', 'CP021', 'CP022', 'CP023', 'CP024', 'CP025', 'CP026', 'CP027', 'CP028', 'CP029', 'CP030', 'CP031', 'CP032', 'CP033', 'CP034', 'CP035', 'CP036', 'CP037', 'CP038', 'CP040', 'CP045', 'CP054', 'CP055', 'CP061', 'CP065', 'CP093', 'CP121', 'CP190']
df = df[df['cp'].isin(cps)]

scenarios = df['scenario'].unique()

fs = 22

for scenario in scenarios:
    df_scenario = df[df['scenario'] == scenario]
    
    # Downsample data thresholds
    df_scenario = df_scenario[df_scenario['threshold'].isin(np.arange(0, 1.1, 0.1))]
    df_scenario = df_scenario.sort_values(by='method')
    
    # Accuracy plot
    plt.figure(figsize=(15, 8))
    sns.boxplot(data=df_scenario, x='threshold', y='accuracy', hue='method')
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
    sns.boxplot(data=df_scenario, x='threshold', y='f1', hue='method')
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
    
    # Calculate mean precision and recall for each method
    for method in df_scenario['method'].unique():
        method_data = df_scenario[df_scenario['method'] == method]
        # Calculate means
        mean_precision = method_data.groupby('threshold')['precision'].mean()
        mean_recall = method_data.groupby('threshold')['accuracy'].mean()
        
        # Create sorted arrays for plotting
        sorted_data = pd.DataFrame({
            'recall': mean_recall,
            'precision': mean_precision
        }).sort_values('recall')
        
        plt.plot(sorted_data['recall'], sorted_data['precision'], 
                marker='o', label=method, linewidth=2, markersize=0)
    
    plt.xlabel('Recall', fontsize=fs)
    plt.ylabel('Precision', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fs, title_fontsize=fs)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{outprefix}_pr_curve_scenario{scenario}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
