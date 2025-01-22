#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
from matplotlib.colors import ListedColormap

def generate_extended_palette(n_colors):
    # Generate a large palette by cycling through seaborn palettes with different hues
    base_palette = sns.color_palette("tab20", min(n_colors, 20))
    extended_palette = base_palette
    while len(extended_palette) < n_colors:
        additional_palette = sns.color_palette("hsv", n_colors)
        extended_palette.extend(additional_palette[:n_colors - len(extended_palette)])
    
    # Remove grey-like colors from the palette
    def is_grey(color):
        r, g, b = color
        return abs(r - g) < 0.1 and abs(g - b) < 0.1 and abs(r - b) < 0.1

    extended_palette = [color for color in extended_palette if not is_grey(color)]
    
    # Ensure the palette has enough colors
    while len(extended_palette) < n_colors:
        additional_palette = sns.color_palette("husl", n_colors)
        extended_palette.extend(additional_palette[:n_colors - len(extended_palette)])
        extended_palette = [color for color in extended_palette if not is_grey(color)]
    
    return extended_palette[:n_colors]

def main(barcode_file, tissue_file, primary_tissue, outfile):
    # Read input TSV files
    barcode_df = pd.read_csv(barcode_file, sep="\t", index_col=0)
    tissue_df = pd.read_csv(tissue_file, sep="\t")

    # Convert tissue_df to have 0 or 1 entries for columns of all unique CSV values from tissues column
    unique_tissues = sorted(set(tissue for tissues in tissue_df["tissues"] for tissue in tissues.split(',')))
    unique_tissues = [primary_tissue] + [tissue for tissue in unique_tissues if tissue != primary_tissue]
    num_tissues = len(unique_tissues)
    start_value = 1000000
    for i, tissue in enumerate(unique_tissues):
        tissue_df[tissue] = tissue_df["tissues"].apply(lambda x: start_value + i if tissue in x.split(',') else 0)

    # Drop the original tissues column
    tissue_df = tissue_df.drop(columns=["tissues"])

    # Make group name the index
    tissue_df = tissue_df.set_index("group_name")

    # Merge the barcode and tissue dataframes
    num_barcode_sites = barcode_df.shape[1]
    # Add an empty buffer column
    barcode_df.insert(num_barcode_sites, '', 0)
    barcode_df = barcode_df.merge(tissue_df, left_index=True, right_index=True)

    # Prepare unique colors for each integer value
    unique_values = np.unique([val for val in barcode_df.iloc[:, :].values.flatten() if val != -1 and val != 0])
    value_to_index = {val: i for i, val in enumerate(sorted(unique_values))}

    # Adjust colors: 0 as white, -1 as grey, others with an extended color palette
    tissue_sums = tissue_df.sum(axis=0)
    non_zero_tissues = tissue_sums[tissue_sums > 0].index.tolist()
    extended_palette = generate_extended_palette(len(unique_values) - len(non_zero_tissues))
    tissue_colors = ["darkgrey", "blue", "red", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "yellow", "black"]
    subset_tissue_colors = [tissue_colors[i] for i in range(num_tissues) if tissue_df[unique_tissues[i]].sum() > 0]
    extended_palette = extended_palette + subset_tissue_colors
    values = sorted(np.unique(barcode_df.values))
    if 0 in values:
        extended_palette = ["white"] + extended_palette
    if -1 in values:
        extended_palette = ["grey"] + extended_palette
    index_to_color = [extended_palette[i] for i, val in enumerate(values)]

    # Create a custom colormap for the dashed line fill
    cmap = ListedColormap(index_to_color)
    cmap.set_bad(color='black', alpha=0.0)  # Set NaN values to be transparent

    # Map the data values to indices
    indexed_data = barcode_df.iloc[:, :].replace(value_to_index)

    # Set up the figure
    plt.figure(figsize=(15, 6))

    # Plot the matrix with discrete colors and borders
    ax = sns.heatmap(indexed_data, 
                     annot=False,  # Remove annotations from the boxes
                     fmt="d", 
                     cmap=cmap, 
                     cbar=False,  # Disable the color bar
                     linewidths=0.5,  # Add box borders
                     linecolor="black")

    # Add a white box over the buffer column
    buffer_col_index = num_barcode_sites
    ax.add_patch(plt.Rectangle((buffer_col_index, -1), 1, barcode_df.shape[0] + 1, fill=True, facecolor='white', edgecolor='none'))

    # Update y-axis labels for the matrix
    plt.yticks(np.arange(len(barcode_df.index)) + 0.5, barcode_df.index, rotation=0)
    plt.xlabel("Barcode sites / Tissues")
    plt.ylabel("Cells")
    plt.title("")

    # Set font size for x and y axis tick labels and titles
    fs=18
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title(ax.get_title(), fontsize=fs)
    ax.set_xlabel(ax.get_xlabel(), fontsize=fs)
    ax.set_ylabel(ax.get_ylabel(), fontsize=fs)

    # Show the plot
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

if __name__ == "__main__":

    barcode_file = sys.argv[1]
    tissue_file = sys.argv[2]
    primary_tissue = sys.argv[3]
    outfile = sys.argv[4]

    # barcode_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3/beam/MMUS1544/CP14/downsampled_char_matrix.tsv"
    # tissue_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3/beam/MMUS1544/CP14/downsampled_tissue_labels.tsv"
    # primary_tissue = "PRL"
    # outfile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3/beam/MMUS1544/CP14/char_matrix_with_tissues.pdf"

    main(barcode_file, tissue_file, primary_tissue, outfile)
