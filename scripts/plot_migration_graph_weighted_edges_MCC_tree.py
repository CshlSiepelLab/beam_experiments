#!/usr/bin/env python3

### This script takes in a MCC tree from TreeAnnotator made from the BEAST2 posterior and then collapses the tree to a migration graph with edges weighted based on node tissue location probabilities. We essentially make a graph of all possible routes and then each route intensity is based on probability of occuring in the tree given node probabilities.

import re, sys
import pandas as pd
import numpy as np
import networkx as nx
import ete3
import matplotlib.pyplot as plt
import matplotlib

def remove_bracket_content(match):
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

# consensus_tree_file = sys.argv[1]

consensus_tree_file = "beast_gundem_2015_2_21_24/A10_sym/tissue_tree_with_trait.tree"

with open(consensus_tree_file, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('tree'):
            tree_info = line

# remove tree name and = from tree info to get only the newick
tree = ''.join(tree_info.split(' ')[3:])

# strip tree string to newick with associated dataframe of annotations
bracket_content_pattern = re.compile(r'\[.*?\]')
annotations = []
newick = re.sub(bracket_content_pattern, remove_bracket_content, tree)
annotations = [re.split(r',(?![^{]*})', x.replace("&", "")) for x in annotations]

annotations = [{key: value.replace("{", "").replace("}", "") for trait in annotation for key, value in [trait.split("=")]} for annotation in annotations]

# label nodes in newick with only leaf labels
node_labeled_newick = label_nodes(newick)
node_labels = get_labels(node_labeled_newick)

# make a dictionary for annotations to node labels
annotations_dict = {}
for node in node_labels:
    annotations_dict[node] = annotations[node_labels.index(node)]

annotations_df = pd.DataFrame.from_dict(annotations_dict, orient='index')
location_probs_df = annotations_df.loc[:, ['location.set', 'location.set.prob']]

# read newick into ete3 Tree
tree = ete3.Tree(node_labeled_newick, format=3)

locations = list(annotations_df['location.set'].values)
uniq_locations = list(set([value for location in locations for value in location.split(",")]))

# make adjacency matrix weighted by location probabilities with source as row index names and recipient as column names
weighted_adjacency_matrix = pd.DataFrame(0, index=uniq_locations, columns=uniq_locations)

for node in tree.traverse():
    if node.is_root():
        continue
    node_name = node.name
    parent_name = node.up.name
    node_locs = list(location_probs_df.loc[node_name, 'location.set'].split(','))
    node_locs_probs = list(location_probs_df.loc[node_name, 'location.set.prob'].split(','))
    parent_locs = list(location_probs_df.loc[parent_name, 'location.set'].split(','))
    parent_locs_probs = list(location_probs_df.loc[parent_name, 'location.set.prob'].split(','))
    for parent_loc in parent_locs:
        parent_loc_index = parent_locs.index(parent_loc)
        parent_loc_prob = float(parent_locs_probs[parent_loc_index])
        for node_loc in node_locs:
            node_loc_index = node_locs.index(node_loc)
            node_loc_prob = float(node_locs_probs[node_loc_index])
            joint_prob = parent_loc_prob * node_loc_prob
            weighted_adjacency_matrix.loc[parent_loc, node_loc] = weighted_adjacency_matrix.loc[parent_loc, node_loc] + joint_prob

# remove diagonal entries to not plot self-migrations
num_tissues = len(weighted_adjacency_matrix)
for i in range(0, num_tissues):
    for j in range(0, num_tissues):
        if i == j:
            weighted_adjacency_matrix.iloc[i,j] = 0

# normalize weighted adjacency matrix so the largest is 1
max_value = weighted_adjacency_matrix.values.max()
weighted_adj_norm = weighted_adjacency_matrix / max_value

# make complete graph with all locations
G = nx.MultiDiGraph()
for loc1 in uniq_locations:
    for loc2 in uniq_locations:
        if loc1 != loc2:
            weight = weighted_adj_norm.loc[loc1, loc2]
            G.add_edge(loc1, loc2, weight=weight)

# Plot the graph
# nx.draw_planar(G, with_labels = True, arrows = True, connectionstyle='arc3, rad = 0.1')



# Extract edge weights
edge_weights = [G.get_edge_data(u,v)[0]['weight'] for u, v in G.edges()]

# Create a colormap based on edge weights
cmap = matplotlib.colormaps['Reds']

# Draw the graph with edge colors
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, connectionstyle='arc3, rad = 0.1', edge_color=edge_colors, edge_cmap=cmap, width=2, font_size=10, font_color='black', font_weight='bold')

# Add a colorbar to show the weight gradient
sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalize)
sm.set_array([])
cbar = plt.colorbar(sm, orientation='vertical')
cbar.set_label('Edge Weight')

plt.show()

