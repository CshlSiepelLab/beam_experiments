#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_bayes_factors(file_path, outfile, bin_width=2, threshold=5):
    # Read the data
    df = pd.read_csv(file_path)
    
    # Get the Bayes factor column
    bf_column = 'bf(gtr-random)'
    
    # Calculate counts above and below threshold
    count_above = sum(df[bf_column] >= threshold)
    count_below = sum(df[bf_column] < threshold)
    
    # Create the figure with specific size
    plt.figure(figsize=(12, 6))
    
    # Calculate bin edges
    min_val = np.floor(df[bf_column].min())
    max_val = np.ceil(df[bf_column].max())
    bins = np.arange(min_val, max_val + bin_width, bin_width)
    
    # Create histogram
    plt.hist(df[bf_column], 
            bins=bins, 
            color='black',
            edgecolor='black',
            linewidth=0.5)
    
    # Add vertical line at threshold
    plt.axvline(x=threshold, color='#666666', linestyle='--', linewidth=2)
    plt.text(threshold + 1, plt.gca().get_ylim()[1]*0.95, 
             f'Threshold = {threshold}', 
             rotation=0, va='top')
    
    # Add count annotations
    plt.text(plt.gca().get_xlim()[1]*0.7, plt.gca().get_ylim()[1]*0.8,
             f'Count ≥ {threshold}: {count_above}',
             color='#4CAF50')
    plt.text(plt.gca().get_xlim()[1]*0.7, plt.gca().get_ylim()[1]*0.7,
             f'Count < {threshold}: {count_below}',
             color='#666666')
    
    # Customize x-axis ticks
    min_tick = np.floor(min_val / 5) * 5
    max_tick = np.ceil(max_val / 5) * 5
    xticks = np.arange(min_tick, max_tick + 5, 5)
    if 0 not in xticks and min_tick < 0 and max_tick > 0:
        xticks = np.sort(np.append(xticks, 0))
    plt.xticks(xticks, rotation=45)
    
    # Labels and title
    plt.xlabel('ln(Bayes Factor) from (GTR - Random)')
    plt.ylabel('Count')
    
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