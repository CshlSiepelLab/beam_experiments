#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from ete3 import Tree
import dendropy
from copy import deepcopy
import networkx as nx

# default colors taken from metient method for consistency in visualizations
DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

def get_migration_counts(tree):
    """
    Calculates the migration counts for each migration edge in the given phylogenetic tree.

    Parameters:
    tree (ete3.Tree): The phylogenetic tree object.

    Returns:
    dict: A dictionary containing the migration counts for each migration edge.
    """
    migration_counts = {}
    for node in tree.traverse():
        if node.is_root():
            continue
        else:
            node_tissue = node.name.split("_")[-1]
            parent_tissue = node.up.name.split("_")[-1]
            if node_tissue == parent_tissue:
                continue
            migration = f"{parent_tissue}_{node_tissue}"
            if migration not in migration_counts:
                migration_counts[migration] = 1
            else:
                migration_counts[migration] += 1
    return migration_counts

def optionally_add_origin_migration(tree, counts, primary):
    """
    Adds the origin branch to the migration counts if the root node location is different from the known origin primary tissue.

    Parameters:
    tree (ete3.Tree): The phylogenetic tree object.
    counts (dict): A dictionary containing migration counts.
    primary (str): The known primary tissue for the origin.

    Returns:
    dict: Updated migration counts including the origin branch if applicable.
    """
    counts_copy = counts
    root_node_location = tree.get_tree_root().name.split("_")[-1]
    if root_node_location != primary:
        migration = f"{primary}_{root_node_location}"
        if migration not in counts_copy:
            counts_copy[migration] = 1
        else:
            counts_copy[migration] += 1
    return counts_copy

def remove_zero_length_nodes(tree):
    """
    Removes zero-length nodes from the given phylogenetic tree.

    Parameters:
    tree (ete3.Tree): The phylogenetic tree object.
    """
    for node in tree.internal_nodes():
        if node.edge_length == 0:
            parent = node.parent_node
            if parent is not None:  # prevents root from being lost
                parent.remove_child(node)
                children = node.child_nodes()
                for child in children:
                    parent.add_child(child)

def dendropy_beast_to_ete_newick_with_strict_locations(tree):
    """
    Converts the given dendropy phylogenetic tree to ete3 Newick format with strict location labels.

    Parameters:
    tree (dendropy.Tree): The dendropy phylogenetic tree object.

    Returns:
    ete3.Tree: The converted ete3 phylogenetic tree object.
    """
    tree_copy = deepcopy(tree)
    i = 0
    for node in tree_copy.preorder_node_iter():
        try:
            prediction = node.taxon.label + "_" + node.annotations.get_value('location')
            node.taxon.label = prediction
        except Exception as e:
            prediction = f"node{i}" + "_" + node.annotations.get_value('location')
            i += 1
        node.label = prediction
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("\'", ""), format=3)
    return ete_tree

def get_posterior_prob_migration_graph(posterior_probs, all_inferred_counts):
    """
    Calculates the probabilistic consensus migration graph based on the posterior probabilities and inferred migration counts.

    Parameters:
    posterior_probs (list): A list of posterior probabilities.
    all_inferred_counts (numpy.ndarray): An array of dictionaries containing the inferred migration counts for each tree.

    Returns:
    dict: A dictionary containing the total counts for each migration edge in the probabilistic consensus migration graph.
    """
    total_counts = {}
    for prob, inferred_counts in zip(posterior_probs, all_inferred_counts):
        for pattern, count in inferred_counts.items():
            for num in range(1, count+1):
                edge = f"{pattern}_{num}"
                if edge not in total_counts:
                    total_counts[edge] = prob
                else:
                    total_counts[edge] += prob
    return total_counts


# inputs
beast_posterior_file = sys.argv[1]
primary_tissue=sys.argv[2]
outdir = sys.argv[3]

# # testing
# beast_posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/yang_2022_real_data_8_22_24/metastabayes/3457_Apc_T4/combined.trees"
# primary_tissue="T"
# outdir = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/yang_2022_real_data_8_22_24/metastabayes/3457_Apc_T4"

# set outfile
outfile = os.path.join(outdir, "beast_posterior_probability_migration_graph.pdf")

# process beast posterior result to get probabilistic consensus migration graph
burnin_percent=0.1
beast_tree_list = dendropy.TreeList()
beast_tree_list.read(path=beast_posterior_file, schema="nexus")
num_beast_trees = len(beast_tree_list)
num_discard = round(num_beast_trees * burnin_percent)
beast_tree_list = beast_tree_list[num_discard:]

posteriors = np.array([float(tree.annotations.get_value('posterior')) for tree in beast_tree_list])
all_inferred_counts = np.empty(len(beast_tree_list), dtype=object)

for i, tree in enumerate(beast_tree_list):
    remove_zero_length_nodes(tree)
    inferred_tree = dendropy_beast_to_ete_newick_with_strict_locations(tree)
    inferred_counts = get_migration_counts(inferred_tree)
    inferred_counts = optionally_add_origin_migration(inferred_tree, inferred_counts, primary_tissue)
    all_inferred_counts[i] = inferred_counts

# fit a gaussian kernel density estimate to the posterior values to get a probability density function
pdf = gaussian_kde(posteriors)
posterior_probs = [pdf(posterior)[0] for posterior in posteriors]
total_posterior_prob = sum(posterior_probs)
posterior_probs = [posterior_prob/total_posterior_prob for posterior_prob in posterior_probs]

# obtain the probabilistic consensus migration graph
graph_dict = get_posterior_prob_migration_graph(posterior_probs, all_inferred_counts)
all_tissues = list(set([value for node in graph_dict.keys() for value in node.split("_")[0:2]]))

# plot the graph
custom_colors = DEFAULT_COLORS
custom_colors = {node: color for node, color in zip(all_tissues, custom_colors[0:len(all_tissues)]) if node != primary_tissue}
custom_colors[primary_tissue] = "black"

G = nx.MultiDiGraph()

for node in all_tissues:
    G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0)

for edge, probability in graph_dict.items():
    if probability > 0.5:   # NEED TO FIX TO REPRESENT THE PROBABILITIES AND NOT A STRICT THRESHOLD
        source, target, num = edge.split('_')
        # # to make the edges transparent based on the probability (this does not work well)
        # source_color = list(mcolors.to_rgba(custom_colors[source]))
        # target_color = list(mcolors.to_rgba(custom_colors[target]))
        # transparency = round(float(probability), 4)
        # source_color[3] = transparency
        # target_color[3] = transparency
        # source_color = "#%02x%02x%02x%02x" % tuple(int(c * 255) for c in source_color)
        # target_color = "#%02x%02x%02x%02x" % tuple(int(c * 255) for c in target_color)
        G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3)


dot = nx.nx_pydot.to_pydot(G)

dot.write_pdf(outfile)


