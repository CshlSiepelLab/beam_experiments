#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_bayes_factors(file_path, outfile, bin_width=1, threshold=5):
    # Read the data
    df = pd.read_csv(file_path)

    # add columns for the rates to the df
    df["mutation_rate"] = [
        name.split("_")[3].replace("mut", "0.") for name in df["name"]
    ]
    df["migration_rate"] = [
        name.split("_")[2].split("/")[1].replace("mig", "1e-") for name in df["name"]
    ]
    unique_mig_rates = sorted(df["migration_rate"].unique(), reverse=True)
    unique_mut_rates = sorted(df["mutation_rate"].unique())
    num_mig_rates = len(unique_mig_rates)
    num_mut_rates = len(unique_mut_rates)

    # Get the Bayes factor column
    bf_column = "bf(gtr-random)"

    # Plot histograms
    plt.figure(figsize=(12, 6))
    fs = 22
    fig, axes = plt.subplots(
        num_mig_rates, num_mut_rates, figsize=(12, 12), sharex=True, sharey=True
    )

    for i, mig in enumerate(unique_mig_rates):
        for j, mut in enumerate(unique_mut_rates):
            subset = df[(df["migration_rate"] == mig) & (df["mutation_rate"] == mut)]
            if subset.empty:
                continue

            ax = axes[i, j]
            ax.hist(subset[bf_column], bins=20, color="grey")
            ax.axvline(x=threshold, color="black", linestyle="--", linewidth=2)
            ax.axvline(x=-threshold, color="black", linestyle="--", linewidth=2)

            ax.set_xlabel("", fontsize=fs)
            ax.set_ylabel("", fontsize=fs)
            ax.tick_params(axis="x", labelsize=fs)
            ax.tick_params(axis="y", labelsize=fs)

            if i == 0:
                ax.set_title(f"Mut rate\n{mut}", fontsize=fs)
            if j == num_mut_rates - 1:
                ax.yaxis.set_label_position("right")
                ax.set_ylabel(
                    f"Mig rate\n{mig}", fontsize=fs, rotation=270, labelpad=50
                )

    fig.supxlabel("ln(Bf)\nfrom ln(ml H_a) - ln(ml H_0)", fontsize=fs)
    fig.supylabel("Count", fontsize=24)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


if __name__ == "__main__":
    file_path = sys.argv[1]
    outfile = sys.argv[2]
    fig = plot_bayes_factors(file_path, outfile)