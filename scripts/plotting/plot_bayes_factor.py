
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_bayes_factors(file_path, outfile, bin_width=1, threshold=1.1):
    # Read the data
    df = pd.read_csv(file_path)

    # Get the Bayes factor column
    bf_column = "bf(gtr-random)"
    # bf_column = 'bf(gtr-gtrNoRLdirectSeeding)'
    # bf_column = 'bf(reseeding-no_reseeding)'

    # Create the figure with specific size
    plt.figure(figsize=(12, 6))

    # Calculate bin edges
    min_val = np.floor(df[bf_column].min())
    max_val = np.ceil(df[bf_column].max())
    bins = np.arange(min_val, max_val + bin_width, bin_width)

    # Create histogram
    plt.hist(df[bf_column], bins=bins, color="grey", edgecolor="grey", linewidth=0)

    # Add vertical line at threshold
    plt.axvline(x=threshold, color="black", linestyle="--", linewidth=2)
    plt.axvline(x=-threshold, color="black", linestyle="--", linewidth=2)

    # Add count annotations
    fs = 20
    count_above_threshold = sum(df[bf_column] > threshold)
    count_below_threshold = sum(df[bf_column] < -threshold)
    count_middle = sum((df[bf_column] > -threshold) & (df[bf_column] < threshold))
    plt.title(
        f"Count > {threshold}: {count_above_threshold}\nCount < {-threshold}: {count_below_threshold}\nCount between {-threshold} and {threshold}: {count_middle}",
        fontsize=fs,
    )

    min_val = -40
    max_val = 120

    increment = 20
    xticks = np.arange(
        np.floor(min_val / increment) * increment,
        np.ceil(max_val / increment) * increment + increment,
        increment,
    )
    if 0 not in xticks:
        xticks = np.sort(np.append(xticks, 0))
    plt.xlim(min_val, max_val)
    plt.xticks(xticks, fontsize=fs)
    plt.xticks(xticks, rotation=0, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlabel("ln(Bf)\nfrom ln(ml H_a) - ln(ml H_0)", fontsize=fs,)
    plt.ylabel("Count", fontsize=fs)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


if __name__ == "__main__":

    file_path = sys.argv[1]
    outfile = sys.argv[2]

    fig = plot_bayes_factors(file_path, outfile)
