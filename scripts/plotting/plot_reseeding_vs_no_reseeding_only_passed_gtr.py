
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_bayes_factors(file_path1, file_path2, outfile, bin_width=1, threshold=3):
    df1 = pd.read_csv(file_path1)
    df2 = pd.read_csv(file_path2)

    # Filter for rows with passed GTR based on threshold
    bf_column1 = "bf(gtr-random)"
    df1 = df1[df1[bf_column1] >= threshold]
    df2 = df2[df2['name'].isin(df1['name'])]

    # Get the Bayes factor column
    bf_column2 = 'bf(reseeding-no_reseeding)'

    # Create the figure with specific size
    plt.figure(figsize=(12, 6))

    # Calculate bin edges
    min_val = np.floor(df2[bf_column2].min())
    max_val = np.ceil(df2[bf_column2].max())
    bins = np.arange(min_val, max_val + bin_width, bin_width)

    # Create histogram
    plt.hist(df2[bf_column2], bins=bins, color="grey", edgecolor="grey", linewidth=0)

    # Add vertical line at threshold
    plt.axvline(x=threshold, color="black", linestyle="--", linewidth=2)
    plt.axvline(x=-threshold, color="black", linestyle="--", linewidth=2)

    # Add count annotations
    fs = 20
    count_above_threshold = sum(df2[bf_column2] > threshold)
    count_below_threshold = sum(df2[bf_column2] < -threshold)
    count_middle = sum((df2[bf_column2] > -threshold) & (df2[bf_column2] < threshold))
    plt.title(
        f"Count > {threshold}: {count_above_threshold}\nCount < {-threshold}: {count_below_threshold}\nCount between {-threshold} and {threshold}: {count_middle}",
        fontsize=fs,
    )

    # min_val = -80
    # max_val = 120
    # plt.xlim(min_val, max_val)

    # Customize x-axis ticks
    increment = 20
    xticks = np.arange(
        np.floor(min_val / increment) * increment,
        np.ceil(max_val / increment) * increment + increment,
        increment,
    )
    if 0 not in xticks:
        xticks = np.sort(np.append(xticks, 0))
    plt.xticks(xticks, rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlabel('ln(Bayes factor)\nfrom ln(Marginal likelihood Reseeding) - ln(Marginal likelihood No Reseeding)', fontsize=fs)
    plt.ylabel("Count", fontsize=fs)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


if __name__ == "__main__":

    file_path1 = sys.argv[1] # gtr vs random
    file_path2 = sys.argv[2] # reseeding vs no reseeding
    outfile = sys.argv[3]

    fig = plot_bayes_factors(file_path1, file_path2, outfile)
