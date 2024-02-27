#!/usr/bin/env python3

### This script reads in a BEAST posterior file for trees labeld with tissues at each node and then subsets all sampled trees for a specified number of trees with the highest posterior probability and then collapses these to migraiton graphs.

import re, sys
import numpy as np
import networkx as nx
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import ete3
import pandas as pd
import matplotlib

def indices_of_top_closest_values(lst, key, n):
    sorted_indices = sorted(range(len(lst)), key=lambda i: abs(lst[i] - key))[:n]
    return sorted_indices

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

def tree_to_migration_graph(tree, primary_tissue, ax, n):
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

    G = nx.MultiDiGraph()

    for node in tree_ete.traverse():
        if node.is_root():
            continue
        node_name = node.name
        parent_name = node.up.name
        node_loc = str(annotations_df.loc[node_name, 'location'])
        parent_loc = str(annotations_df.loc[parent_name, 'location'])
        if parent_loc != node_loc:
            G.add_edge(parent_loc, node_loc)

    # Draw the graph with edge colors
    nodes = sorted(list(G.nodes()))
    G = G.subgraph(nodes)
    node_colors = range(len(nodes))
    node_cmap = matplotlib.cm.get_cmap('tab20', len(nodes))

    # Find the node corresponding to the primary tissue
    primary_tissue_node = [node for node in G.nodes() if node == primary_tissue][0]

    # Create positions for the nodes
    max_width = ax.get_position().width
    pos = {}
    row_height = 0.1
    num_nodes = len(nodes)

    for i, node in enumerate(G.nodes()):
        if node == primary_tissue:
            pos[node] = (max_width / 2, 0)
        else:
            pos[node] = ((max_width / num_nodes) * (i + 0.5), -row_height) 

    nx.draw(G, 
            pos=pos, 
            ax=ax, 
            with_labels=False,
            font_size = 10, 
            connectionstyle='arc3, rad = 0.2', 
            arrowsize = 20,
            width = 1,
            font_color='black', 
            font_weight='bold', 
            node_shape = 's',
            node_size = 1000,
            node_color = node_colors,
            cmap=node_cmap)

    legend_labels = {loc: node_cmap(i) for i, loc in enumerate(list(G.nodes()))}
    legend_handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color, markersize=10) for color in legend_labels.values()]
    ax.legend(legend_handles, legend_labels.keys(), title='Node locations', loc='upper left', bbox_to_anchor=(0.7, 1))
    plt.tight_layout()




posterior_file = sys.argv[1]
primary_tissue = sys.argv[2]

# posterior_file = "beast_gundem_2015_2_21_24/A10_sym/tissue_tree_with_trait.trees"
# primary_tissue = "prostate"

# set the number of trees to obtain as graphs
n = 3

names_dict = {}
trees = []

with open(posterior_file, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('tree'):
            trees.append(line)
        # lines that begin with a number are translate lines
        elif line and line[0].isdigit():
            key_value = line.split(' ')
            key = key_value[0]
            # remove trailing comma for translate values
            value = key_value[1][0:-1]
            names_dict[key] = value

# sort trees by posterior
pattern = re.compile(r'tree STATE_\d+ = \[&posterior=(-?\d+\.\d+)\]')
sorted_trees = sorted(trees, key=lambda s: float(pattern.search(s).group(1)), reverse = True)

# get the top n trees with highest posterior
top_n_trees = sorted_trees[0:n]

# get the top n maximum clade credibility trees by finding closest to the peak of probability density function
posterior_values = re.findall(r'\[&posterior=(-?\d+\.\d+)\]', "".join(trees))
posterior_values = [round(float(value), 2) for value in posterior_values]

kde = gaussian_kde(posterior_values)
x_values = np.linspace(min(posterior_values), max(posterior_values), 1000)
peak_value = x_values[np.argmax(kde(x_values))]
peak_density = kde(peak_value)

# # plot posterior values to see peak
# plt.plot(x_values, kde(x_values), label='Posterior Density')
# plt.hist(posterior_values, bins=100, density=True, alpha=0.5, color='green', label='Posterior Histogram')
# plt.scatter(peak_value, peak_density, color='red', label=f'Peak: {peak_value:.2f}')
# plt.xlabel("Posterior")
# plt.ylabel("Density")
# plt.show()

# get n number of trees closest to the peak value of the posterior density
top_n_indices = indices_of_top_closest_values(posterior_values, peak_value, n)
mcc_top_n_trees = [trees[i] for i in top_n_indices]



# for tree in top_n_trees:
fig, axes = plt.subplots(2, n, figsize=(12, 6))

i = 0
for tree in top_n_trees:
    tree_to_migration_graph(tree, primary_tissue, axes[0][i], i)
    i = i + 1

i = 0
for tree in mcc_top_n_trees:
    tree_to_migration_graph(tree, primary_tissue, axes[1][i], i)
    i = i + 1

cols = list(range(1, n+1, 1))
rows = ["Highest posterior", "Highest density"]
for i, axis in enumerate(axes):
    prefix = rows[i]
    for ax, col in zip(axis, cols):
        name = prefix + " " + str(col)
        ax.set_title(name, fontsize=16)

plt.tight_layout()
# plt.show()

outfile = posterior_file.split(".")[0] + "_migration_graphs.pdf"
plt.savefig(outfile)


