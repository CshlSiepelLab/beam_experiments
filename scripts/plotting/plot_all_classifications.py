
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def preprocess_data(data):
    # Group by threshold and calculate the percentage of True values for each category
    grouped = (
        data.groupby("threshold")
        .agg(
            met_to_met=("met_to_met", lambda x: (x.sum() / len(x)) * 100),
            met_to_primary=("met_to_primary", lambda x: (x.sum() / len(x)) * 100),
        )
        .reset_index()
    )

    # Sort data by threshold in ascending order
    grouped = grouped.sort_values(by="threshold", ascending=True)
    return grouped


def plot_data(data, outfile):
    # Plot the data
    fig, ax = plt.subplots(figsize=(12, 6))
    fs = 22

    # Define the width of the bars
    bar_width = 0.35

    # Define the positions of the bars
    index = np.arange(len(data))

    # Plot the bars
    bar1 = ax.bar(
        index,
        data["met_to_met"],
        bar_width,
        # label="BEAM Met to Met",
        label="MACH2 Met to Met",
        color="#4c72b0",
        alpha=1.0,
    )
    bar2 = ax.bar(
        index + bar_width,
        data["met_to_primary"],
        bar_width,
        # label="BEAM Primary Reseeding",
        label="MACH2 Primary Reseeding",
        color="#dd8452",
        alpha=1.0,
    )

    # Set y-axis limits
    ax.set_ylim(0, 100)

    # Add labels and title
    ax.set_xlabel("Posterior probability", fontsize=fs, labelpad=10)
    ax.set_ylabel("Percentage of data (%)", fontsize=fs, labelpad=10)
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(data["threshold"].astype(str), fontsize=fs)
    ax.tick_params(axis="y", labelsize=fs)

    # # Optional: drow a horizontal line at specified y-value
    # ax.axhline(y=54.0669856459, color="#3a5f8a", linestyle="--", linewidth=2, label="MACH2 M2M")    # mach2 m2m average (value obtained seperately)
    # ax.axhline(y=27.7853725222146, color="#b36238", linestyle="--", linewidth=2, label="MACH2 PR")    # mach2 pr average (value obtained seperately)

    # Add legend with improved placement
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 0.75), frameon=False, fontsize=16)

    # Improve layout for publication
    plt.tight_layout(pad=2.0)
    plt.savefig(
        outfile, bbox_inches="tight", dpi=300
    )  # Higher DPI for publication quality
    plt.close()

    # print average across threholds for each met to met and met to primary
    print(data["met_to_met"].mean())
    print(data["met_to_primary"].mean())


def main():
    file_path = sys.argv[1]
    outfile = sys.argv[2]

    data = pd.read_csv(file_path)
    data = preprocess_data(data)
    plot_data(data, outfile)


if __name__ == "__main__":
    main()
