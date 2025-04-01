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

    # Get the Bayes factor column
    bf_column = "bf(gtr-random)"

    # Create the figure with specific size
    plt.figure(figsize=(12, 6))

    # Create histogram
    unique_mig_rates = df["migration_rate"].unique()
    unique_mut_rates = df["mutation_rate"].unique()

    # Determine the number of unique migration and mutation rates
    num_mig_rates = len(unique_mig_rates)
    num_mut_rates = len(unique_mut_rates)

    # Create a figure with subplots
    fig, axes = plt.subplots(
        num_mig_rates, num_mut_rates, figsize=(12, 12), sharex=True, sharey=True
    )

    fs = 22

    # Sort unique rates
    unique_mig_rates = sorted(unique_mig_rates, reverse=True)
    unique_mut_rates = sorted(unique_mut_rates)

    # Plot each subplot
    for i, mig in enumerate(unique_mig_rates):
        for j, mut in enumerate(unique_mut_rates):
            subset = df[(df["migration_rate"] == mig) & (df["mutation_rate"] == mut)]
            if subset.empty:
                continue

            ax = axes[i, j]
            ax.hist(subset[bf_column], bins=20, color="grey")
            ax.axvline(x=threshold, color="black", linestyle="--", linewidth=2)

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

    # Set common labels
    fig.supxlabel("ln(Bayes Factor) from (GTR - Random)", fontsize=fs)
    fig.supylabel("Count", fontsize=24)

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


if __name__ == "__main__":

    file_path = sys.argv[1]
    outfile = sys.argv[2]

    # file_path = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_2_25_25_data_from_8_19_24/marginal_likelihoods.csv"
    # outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_2_25_25_data_from_8_19_24/marginal_likelihoods_variable_rates.pdf"

    fig = plot_bayes_factors(file_path, outfile)
