#!/usr/bin/env python3

import re, sys, os
import random
import pandas as pd
import numpy as np
import networkx as nx
from ete3 import Tree
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import to_pydot
import matplotlib
from IPython.display import Image, display
import matplotlib.image as mpimg
from io import BytesIO

tree_file=sys.argv[1]
# tree_file="/grid/siepel/home_norepl/staklins/stephen_data/beast_bayesian_migration_graph_inference/variable_migration_and_mutation_rates_8_19_24_data_from_8_19_24/raw_data/mig4_mut001_231/tissue_labeled_tree.nwk"

outfile = tree_file.replace(".nwk", "_migration_graph.pdf")

primary_tissue="P"

tree = Tree(tree_file, format=8)

migration_counts = {}

if "_" not in tree.get_tree_root().name:
    tree.get_tree_root().name = f'0_{primary_tissue}'
else:
    # handle the origin branch between the origin and root if the root has a tissue label already
    root_tissue = tree.get_tree_root().name.split("_")[1]
    if root_tissue != primary_tissue:
        migration = f"{root_tissue}_{primary_tissue}"
        migration_counts[migration] = 1

all_tissues = [primary_tissue]

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

        if node_tissue not in all_tissues:
            all_tissues.append(node_tissue)
        if parent_tissue not in all_tissues:
            all_tissues.append(parent_tissue)

all_tissues = sorted(all_tissues)

num_nodes = len(all_tissues)

# default colors taken from metient method for consistency in visualizations
DEFAULT_COLORS = ["#6aa84f", "#be5742e1", "#6fa8dc", "#e69138", "#9e9e9e", "#c27ba0","brown", "black", "darkgreen", "purple", "blue"]*3
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

