#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_bayes_factors(file_path, outfile, bin_width=1, threshold=5):
    # Read the data
    df = pd.read_csv(file_path)
    
    # Get the Bayes factor column
    bf_column = 'bf(gtr-random)'
    
    # Create the figure with specific size
    plt.figure(figsize=(12, 6))
    
    # Calculate bin edges
    min_val = np.floor(df[bf_column].min())
    max_val = np.ceil(df[bf_column].max())
    bins = np.arange(min_val, max_val + bin_width, bin_width)
    
    # Create histogram
    plt.hist(df[bf_column], 
            bins=bins, 
            color='grey',
            edgecolor='grey',
            linewidth=0)
    
    # Add vertical line at threshold
    plt.axvline(x=threshold, color='black', linestyle='--', linewidth=2)
    
    # Add count annotations
    fs = 20
    plt.text(plt.gca().get_xlim()[1]*0.6, 
             plt.gca().get_ylim()[1]*0.8,
             f'Count ≥ {threshold}: {sum(df[bf_column] >= threshold)}',
             color='black', fontsize=fs)
    
    plt.text(plt.gca().get_xlim()[1]*0.6, 
             plt.gca().get_ylim()[1]*0.7,
             f'Count < {threshold}: {sum(df[bf_column] < threshold)}',
             color='black', fontsize=fs)
    
    # Customize x-axis ticks
    increment = 20
    xticks = np.arange(np.floor(min_val / increment) * increment, np.ceil(max_val / increment) * increment + increment, increment)
    if 0 not in xticks:
        xticks = np.sort(np.append(xticks, 0))
    plt.xticks(xticks, rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlabel('ln(Bayes factor)\nfrom ln(Marginal likelihood GTR) - ln(Marginal likelihood Random)', fontsize=fs)
    plt.ylabel('Count', fontsize=fs)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

if __name__ == "__main__":
    
    file_path = sys.argv[1]
    outfile = sys.argv[2]

    # file_path = "marginal_likelihoods.csv"
    # outfile = "bayes_factor_distribution.pdf"

    fig = plot_bayes_factors(file_path, outfile)