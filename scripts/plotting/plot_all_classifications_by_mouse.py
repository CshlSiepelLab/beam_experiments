
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def preprocess_data(data):
    # Add mouse column by splitting 'name' on '_' and taking field 1
    data['mouse'] = data['name'].str.split('_').str[0]
    return data

def summarize_by_mouse_and_threshold(data):
    # Group by mouse and threshold, calculate percentage of True for each category
    grouped = (
        data.groupby(['mouse', 'threshold'])
        .agg(
            met_to_met=("met_to_met", lambda x: (x.sum() / len(x)) * 100),
            met_to_primary=("met_to_primary", lambda x: (x.sum() / len(x)) * 100),
        )
        .reset_index()
    )
    # Sort by mouse and threshold
    grouped = grouped.sort_values(by=['mouse', 'threshold'])
    return grouped

def plot_data(data, outfile):
    mice = sorted(data['mouse'].unique())
    n_mice = len(mice)
    ncols = 4
    nrows = 3
    fs = 18

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4), squeeze=False)
    bar_width = 0.35

    # Set mach2 lines to be drawn by mouse
    mach2_met = {
        "MMUS1457": 34.883721,
        "MMUS1466": 42.307692,
        "MMUS1467": 42.857143,
        "MMUS1469": 34.166667,
        "MMUS1492": 0.000000,
        "MMUS1495": 61.403509,
        "MMUS1544": 69.767442,
        "MMUS1588": 39.130435,
        "MMUS1874": 5.882353,
        "MMUS1875": 9.523810,
    }
    mach2_pr = {
        "MMUS1457": 18.604651,
        "MMUS1466": 7.692308,
        "MMUS1467": 22.857143,
        "MMUS1469": 30.000000,
        "MMUS1492": 28.571429,
        "MMUS1495": 34.085213,
        "MMUS1544": 20.930233,
        "MMUS1588": 17.391304,
        "MMUS1874": 5.882353,
        "MMUS1875": 14.285714,
    }

    for idx, mouse in enumerate(mice):
        ax = axes[idx // ncols, idx % ncols]
        mouse_data = data[data['mouse'] == mouse]
        index = np.arange(len(mouse_data))

        bar1 = ax.bar(
            index,
            mouse_data["met_to_met"],
            bar_width,
            label="BEAM Met to Met",
            color="#4c72b0",
            alpha=1.0,
        )
        bar2 = ax.bar(
            index + bar_width,
            mouse_data["met_to_primary"],
            bar_width,
            label="BEAM Primary Reseeding",
            color="#dd8452",
            alpha=1.0,
        )

        ax.set_ylim(0, 100)
        ax.set_xlabel("Threshold", fontsize=fs)
        ax.set_ylabel("Percent (%)", fontsize=fs)
        ax.set_xticks(index + bar_width / 2)
        ax.set_xticklabels(mouse_data["threshold"].astype(str), fontsize=fs-2)
        ax.tick_params(axis="y", labelsize=fs-2)
        ax.set_title(f"Mouse {mouse}", fontsize=fs+2)
        if idx == 0:
            ax.legend(fontsize=fs-2)

        # Add mach2 lines
        ax.axhline(y=mach2_met[mouse], color="#3a5f8a", linestyle="--", linewidth=2, label="MACH2 Met to Met")    # mach2 met to met average (value obtained seperately)
        ax.axhline(y=mach2_pr[mouse], color="#b36238", linestyle="--", linewidth=2, label="MACH2 Primary Reseeding")    # mach2 primary reseeding average (value obtained seperately)


    # Hide unused subplots
    for idx in range(n_mice, nrows * ncols):
        fig.delaxes(axes[idx // ncols, idx % ncols])

    plt.tight_layout(pad=2.0)
    plt.savefig(outfile, bbox_inches="tight", dpi=300)
    plt.close()

def main():
    file_path = sys.argv[1]
    outfile = sys.argv[2]

    data = pd.read_csv(file_path)
    data = preprocess_data(data)
    data = summarize_by_mouse_and_threshold(data)
    plot_data(data, outfile)

if __name__ == "__main__":
    main()
