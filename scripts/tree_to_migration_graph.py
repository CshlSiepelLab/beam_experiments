#!/usr/bin/env python3

import re, sys, os
import random
import pandas as pd
import numpy as np
import networkx as nx
from ete3 import Tree
import matplotlib.pyplot as plt
import matplotlib

# tree_file=sys.argv[1]
tree_file="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/machina_tree_all_tissue_labels.nwk"
primary_tissue="P"

tree = Tree(tree_file, format=8)

tree.get_tree_root().name = f'0_{primary_tissue}'

migration_counts = {}

for node in tree.traverse():
    if node.is_root() or node.is_leaf():
        continue
    else:
        node_tissue = node.name.split("_")[1]
        parent_tissue = node.up.name.split("_")[1]
        if node_tissue == parent_tissue:
            continue
        migration = f"{parent_tissue}_{node_tissue}"
        if migration not in migration_counts:
            migration_counts[migration] = 1
        else:
            migration_counts[migration] += 1

G = nx.MultiDiGraph()

for edge, count in migration_counts.items():
    source, target = edge.split('_')
    for _ in range(count):
        G.add_edge(source, target)


# plot graph
fig, ax = plt.subplots(figsize=(8, 8))
max_width = ax.get_position().width
pos = {}
row_height = 0.05
num_nodes = len(G.nodes())

for i, node in enumerate(G.nodes()):
    if node == primary_tissue:
        pos[node] = (max_width / 2, 0)
    else:
        pos[node] = ((max_width / num_nodes) * (i + 0.5), -row_height + random.uniform(0, 0.025)) 


node_colors = ["black", "red", "green", "blue", "orange", "purple", "brown", "pink", "gray", "gold"]
node_colors = node_colors[0:num_nodes]

nx.draw(G, 
        pos=pos, 
        ax=ax, 
        with_labels=False, 
        connectionstyle='arc3, rad = 0.2', 
        arrowsize = 20,
        font_size=10, 
        font_color='black', 
        font_weight='bold', 
        node_shape = 's',
        node_size = 1000,
        node_color = node_colors)

legend_labels = {loc: node_color for loc, node_color in zip(G.nodes(), node_colors)}
legend_handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color, markersize=10) for color in legend_labels.values()]
ax.legend(legend_handles, legend_labels.keys(), title='Node Locations', loc='upper left', bbox_to_anchor=(0.9, 1))

# save graph to file
outfile = tree_file.split(".")[0] + "_migration_graph.pdf"
plt.savefig(outfile)

#plt.show()