
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
        "MMUS1457": 70.83333333333336,
        "MMUS1466": 70.58823529411767,
        "MMUS1467": 57.692307692307686,
        "MMUS1469": 37.61467889908257,
        "MMUS1492": 0.0,
        "MMUS1495": 65.38461538461537,
        "MMUS1544": 83.09859154929576,
        "MMUS1588": 47.36842105263157,
        "MMUS1874": 13.043478260869565,
        "MMUS1875": 12.5,
    }
    mach2_pr = {
        "MMUS1457": 33.33333333333332,
        "MMUS1466": 11.76470588235294,
        "MMUS1467": 30.76923076923077,
        "MMUS1469": 28.440366972477072,
        "MMUS1492": 44.44444444444445,
        "MMUS1495": 37.5,
        "MMUS1544": 24.144869215291745,
        "MMUS1588": 10.526315789473683,
        "MMUS1874": 8.695652173913045,
        "MMUS1875": 18.75,
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
            ax.legend(fontsize=fs-12)

        # Add mach2 lines
        ax.axhline(y=mach2_met[mouse], color="#3a5f8a", linestyle="--", linewidth=2, label="MACH2 Met to Met")    # mach2 met to met average (value obtained seperately)
        ax.axhline(y=mach2_pr[mouse], color="#b36238", linestyle="--", linewidth=2, label="MACH2 Primary Reseeding")    # mach2 primary reseeding average (value obtained seperately)
        
        print(mouse)
        print("m2m", mouse_data["met_to_met"].mean())
        print("pr", mouse_data["met_to_primary"].mean())


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
