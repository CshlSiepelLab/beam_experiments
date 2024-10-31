#!/usr/bin/env python3

import sys
import pandas as pd
from ete3 import Tree
from itertools import combinations
import matplotlib.pyplot as plt
import random
from ete3 import Tree
import plot_phylo


# inputs
char_matrix_file = sys.argv[1]
tree_file = sys.argv[2]
tissue_labels_file = sys.argv[3]
thresh = float(sys.argv[4])
outprefix = sys.argv[5]

# char_matrix_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/raw_data/mS_854/mS_854_indel_character_matrix.tsv"
# tree_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/laml/mS_854/mS_854_laml_trees.nwk"
# tissue_labels_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/raw_data/mS_854/cell_tree_seed1833437564.labeling"
# thresh = 1.5 # to be adjusted as a scalar of the tree height
# outprefix = "test"

# whether to make plots to show the downsampling result compared to the original tree
plot = True

# read in tissue data
tissue_labels_original = pd.read_csv(tissue_labels_file, sep=r'\s+', index_col=None, names=['group_name', 'tissues'], dtype=str)
tissue_group_names_original = [str(l) for l in tissue_labels_original['group_name'].to_list()]

# only keep clones with one tissue label
tissue_labels = tissue_labels_original[tissue_labels_original['tissues'].str.contains(',', na=False) == False]
tissue_group_names = [str(l) for l in tissue_labels['group_name'].to_list()]

# read in character matrix
char_matrix = pd.read_csv(char_matrix_file, sep='\t', index_col=0)
char_matrix.index = char_matrix.index.astype(str)

# read in tree
try:
    tree = Tree(tree_file, format=3)
except:
    tree = Tree(tree_file, format=5)

# get the actual threshold real value as a distance based on the tree height, assuming the tree is ultrametric
farthest_leaf, total_length = tree.get_farthest_leaf()
threshold = thresh * total_length * 2
print("Real branch length threshold value: ", threshold)

# get leaf nodes and their names
tree_leaf_nodes = []
for node in tree.traverse():
    if node.is_leaf():
        tree_leaf_nodes.append(node)

tree_leaves = [node.name for node in tree_leaf_nodes]

# run basic checks
assert set(tree_leaves) == set(tissue_group_names_original), "The tree leaves from the newick file input do not match the tissue group names input."

num_clones = len(tree_leaves)
print("Initial number of clones: ", num_clones)
print("Number of clones with more than one tissue label: ", len(tissue_group_names_original) - len(tissue_group_names))
print("Number of clones with one tissue label: ", len(tissue_group_names))

# downsample the tree tips by considering pairwise distances between candidate tips with one tissue label while leaving tips with more than one tissue label untouched
candidates = [node for node in tree_leaf_nodes if str(node.name) in tissue_group_names]

nodes_to_remove = set()
distance_values = []
for tissue in tissue_labels['tissues'].unique():
    nodes_in_tissue = [node for node in candidates if tissue_labels.loc[tissue_labels['group_name'] == str(node.name), 'tissues'].values[0] == str(tissue)]
    distances=[[node1, node2, node1.get_distance(str(node2.name))] for node1, node2 in combinations(nodes_in_tissue, 2)]
    distance_values.extend(distances)
    for node1, node2, distance in distances:
        if distance < threshold:
            if node1 not in nodes_to_remove and node2 not in nodes_to_remove:
                char_matrix_node1 = char_matrix.loc[node1.name]
                char_matrix_node2 = char_matrix.loc[node2.name]
                set1 = set(char_matrix_node1[char_matrix_node1 != 0][char_matrix_node1 != -1])
                set2 = set(char_matrix_node2[char_matrix_node2 != 0][char_matrix_node2 != -1])
                
                # if they have the same number of unique mutations, choose the one with the least missing data
                if len(set1) == len(set2):
                    node1_missing = (char_matrix_node1 == -1).sum()
                    node2_missing = (char_matrix_node2 == -1).sum()
                    if node1_missing > node2_missing:
                        nodes_to_remove.add(node1)
                    elif node2_missing > node1_missing:
                        nodes_to_remove.add(node2)
                    # if they have the same number of unique mutations and missing data, choose randomly
                    else:
                        nodes_to_remove.add(random.choice([node1, node2]))
                # main priority is to keep the clone with the most unique mutations
                elif len(set1) > len(set2):
                    nodes_to_remove.add(node2)
                else:
                    nodes_to_remove.add(node1)
                
# Plot histogram of values and threshold
if plot:
    plt.hist([distance for _, _, distance in distance_values], bins=100, edgecolor='black', color='grey')
    plt.axvline(x=threshold, color='r', linestyle='dashed', linewidth=1)
    plt.title('Histogram of pairwise distances for candidate tips')
    plt.xlabel('Distance (branch length)')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(outprefix + "_downsample_threshold_histogram.png")
    plt.close()

# remove nodes from the tree in a copy of the tree
tree_copy = tree.copy()
remove_names = [node.name for node in nodes_to_remove]
keep_names = [name for name in tree_leaves if name not in remove_names]
tree_copy.prune(keep_names, preserve_branch_length=True)

# filter tissue labels
tissue_labels_filtered = tissue_labels_original[~tissue_labels_original['group_name'].isin(remove_names)]

# filter char matrix
char_matrix_filtered = char_matrix.drop(remove_names)

# assign colors to the tip names based on the tissue labels
DEFAULT_COLORS = ["black", "red", "green", "purple", "blue", "orange", "brown", "pink", "grey", "yellow", "cyan"]
tissue_colors = {}
i=0
for tissue in tissue_labels['tissues'].unique():
    tissue_colors[tissue] = DEFAULT_COLORS[i]
    i += 1

color_dict = {}
for node in tree.traverse():
    if node.is_leaf():
        if str(node.name) in tissue_group_names:
            color_dict[str(node.name)] = tissue_colors[tissue_labels.loc[tissue_labels['group_name'] == str(node.name), 'tissues'].values[0]]

# plot the trees side by side
if plot:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    plot_phylo.plot_phylo(tree.write(), ax=ax1, col_dict=color_dict, show_support=False)
    ax1.set_title(f"Original Tree: {num_clones} tips")
    plot_phylo.plot_phylo(tree_copy.write(), ax=ax2, col_dict=color_dict, show_support=False, reverse=True) 
    ax2.set_title(f"Downsampled Tree: {len(keep_names)} tips")
    ax2.axvline(x=(ax2.get_xlim()[1] * thresh), color='r', linestyle='dashed', linewidth=1)
    fig.suptitle(f"Threshold at {thresh*100}% of tree height")
    plt.tight_layout()
    plt.savefig(outprefix + "_downsampled_tree_comparison.png")
    plt.close()

# output downsampled tree and tissue labels
tree_copy.write(outfile=f"{outprefix}_downsampled_tree.nwk", format=5)

# tissue_labels_filtered.to_csv(outprefix + "_tissue_labels.tsv", sep=',', header=False)  # for sim data
tissue_labels_filtered.to_csv(outprefix + "_tissue_labels.tsv", sep='\t', header=True, index=False) # for yang data

# output downsampled character matrix
char_matrix_filtered.to_csv(outprefix + "_char_matrix.tsv", sep='\t')
