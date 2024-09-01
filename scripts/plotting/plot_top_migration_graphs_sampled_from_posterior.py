#!/usr/bin/env python3

import re
import sys
import os
import numpy as np
import networkx as nx
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import ete3
import pandas as pd
import matplotlib
import random
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import to_pydot

# default colors taken from metient method for consistency in visualizations
DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3

def remove_bracket_content(match):
    global annotations
    annotations.append(match.group()[1:-1])
    return ''

def get_labels(newick):
    label_pattern = re.compile(r'([A-Za-z0-9_]+):')
    leaf_labels = label_pattern.findall(newick)
    return leaf_labels

def label_nodes(newick):
    leaf_labels = get_labels(newick)
    try:
        leaf_labels_max = max(map(int, leaf_labels))
    except ValueError:
        # Handle the case where leaf labels are not convertible to integers
        leaf_labels_max = len(leaf_labels)
    # add labels to nodes
    start_label = leaf_labels_max + 1
    node_labeled_newick = ""
    parts = newick.split(")")
    for part in parts[:-1]:
        part = part + ")" + str(start_label)
        node_labeled_newick += part
        start_label += 1
    node_labeled_newick += parts[-1]
    return node_labeled_newick

def tree_to_migration_graph(tree, primary_tissue, outfile):
    global annotations
    tree = ''.join(tree.split(' ')[4:])
    bracket_content_pattern = re.compile(r'\[.*?\]')
    annotations = []
    newick = re.sub(bracket_content_pattern, remove_bracket_content, tree)
    node_labeled_newick = label_nodes(newick)
    node_labels = get_labels(node_labeled_newick)
    annotations = [re.split(r',(?![^{]*})', x.replace("&", "")) for x in annotations]
    annotations = [{key: value.replace("{", "").replace("}", "") for trait in annotation for key, value in [trait.split("=")]} for annotation in annotations]
    annotations_dict = {}
    for node in node_labels:
        annotations_dict[node] = annotations[node_labels.index(node)]
    annotations_df = pd.DataFrame.from_dict(annotations_dict, orient='index')
    annotations_df.replace('"', '', regex=True, inplace=True)

    tree_ete = ete3.Tree(node_labeled_newick, format=3)

    migration_counts = {}
    all_tissues = [primary_tissue]

    for node in tree_ete.traverse():
        if node.is_leaf():
            continue
        elif node.is_root():
            node_name = node.name
            node_loc = str(annotations_df.loc[node_name, 'location'])
            if node_loc != primary_tissue:
                migration = f"{node_loc}_{primary_tissue}"
                migration_counts[migration] = 1
        else:
            node_tissue = str(annotations_df.loc[node.name, 'location'])
            parent_tissue = str(annotations_df.loc[node.up.name, 'location'])
            if node_tissue == parent_tissue:
                continue
            migration = f"{parent_tissue}_{node_tissue}"
            if migration not in migration_counts:
                migration_counts[migration] = 1
            else:
                migration_counts[migration] += 1

            if node_tissue not in all_tissues:
                all_tissues.append(node_tissue)
            if parent_tissue not in all_tissues:
                all_tissues.append(parent_tissue)

        all_tissues = sorted(all_tissues)
        num_nodes = len(all_tissues)

        custom_colors = DEFAULT_COLORS

        custom_colors = {node: color for node, color in zip(all_tissues, custom_colors[0:num_nodes]) if node != primary_tissue}
        custom_colors[primary_tissue] = "black"

        G = nx.MultiDiGraph()

        for node in all_tissues:
            G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0)

        for edge, count in migration_counts.items():
            source, target = edge.split('_')
            for _ in range(count):
                G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3)

        dot = nx.nx_pydot.to_pydot(G)
        dot.write_pdf(outfile)


# inputs
posterior_file = sys.argv[1]
primary_tissue = sys.argv[2]
n = int(sys.argv[3])

# # testing
# posterior_file = "/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/variable_migration_and_mutation_rates_8_19_24_data_from_8_19_24/metastabayes/mig6_mut005_16247/combined.trees"
# primary_tissue = "P"
# n = 3


names_dict = {}
trees = []

with open(posterior_file, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('tree'):
            trees.append(line)
        # lines that begin with a number and have two fields are translate lines
        elif line and line[0].isdigit() and len(line.split(' ')) > 1:
            key_value = line.split(' ')
            key = key_value[0]
            # remove trailing comma for translate values
            value = key_value[1].replace(',', '')
            names_dict[key] = value

# sort trees by posterior
pattern = re.compile(r'tree STATE_\d+ = \[&posterior=(-?\d+\.\d+),')
sorted_trees = sorted(trees, key=lambda s: float(pattern.search(s).group(1)), reverse = True)
posterior_values = re.findall(r'\[&posterior=(-?\d+\.\d+),', "".join(sorted_trees))
posterior_values = [round(float(value), 2) for value in posterior_values]

# get peak value and density
kde = gaussian_kde(posterior_values)
x_values = np.linspace(min(posterior_values), max(posterior_values), 1000)
max_x = max(x_values)
max_x_y = kde(max_x)
peak_value = x_values[np.argmax(kde(x_values))]
peak_density = kde(peak_value)

# # plot posterior values to see peak
# plt.plot(x_values, kde(x_values))
# plt.hist(posterior_values, bins=100, density=True, alpha=0.5, color='green')
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.xlabel("Posterior", fontsize=24)
# plt.ylabel("Density", fontsize=24)
# plt.show()
# plt.close()

# get n number of trees closest to the peak value of the posterior density
mcc_top_n_trees = [trees[i] for i in sorted(range(len(posterior_values)), key=lambda i: abs(posterior_values[i] - peak_value))[:n]]

# sample trees, convert to migration graphs, and output as files
i = 1
for tree in mcc_top_n_trees:
    outfile = posterior_file.split(".")[0] + "_migration_graph_" + str(i) + ".pdf"
    tree_to_migration_graph(tree, primary_tissue, outfile)
    i = i + 1
