#!/usr/bin/env python3

import re
import sys
import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from ete3 import Tree
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
    newick_str = ''.join(tree.split(' ')[4:])
    tree = Tree(newick_str, format=1)
    all_tissues = set()
    all_tissues.add(primary_tissue)
    i = 1
    for node in tree.traverse():
        # split the actual name from the annotation read in from beast where the node.name has the form of 7[&location="TBL"] for a tip or [&location="TBL"] for a internal node without a name
        node.name, annotation = node.name.split("[")
        node.tissue = re.search(r'&location="([^"]+)"', annotation).group(1)

        # Label internal nodes without a name
        if node.name == "":
            if node.is_root():
                node.name = "root"
            else:
                node.name = f"node{i}"
                i += 1

        if node.tissue not in all_tissues:
            all_tissues.add(node.tissue)

    migration_counts = {}

    for node in tree.traverse():
        if node.is_root() and node.tissue != primary_tissue:
                migration = f"{primary_tissue}_{node.tissue}"
                if migration in migration_counts:
                    migration_counts[migration] += 1
                else:
                    migration_counts[migration] = 1
        elif not node.is_root() and node.up.tissue != node.tissue:
                migration = f"{node.up.tissue}_{node.tissue}"
                if migration in migration_counts:
                    migration_counts[migration] += 1
                else:
                    migration_counts[migration] = 1
        
    all_tissues = [primary_tissue] + sorted(list(all_tissues - {primary_tissue}))
    num_nodes = len(all_tissues)

    custom_colors = DEFAULT_COLORS

    custom_colors = {node: color for node, color in zip(all_tissues, custom_colors[0:num_nodes]) if node != primary_tissue}
    custom_colors[primary_tissue] = "black"

    G = nx.MultiDiGraph()

    for node in all_tissues:
        G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0)

    for edge, count in migration_counts.items():
        if count == 1:
            label = ""
        else:
            label = str(count)
        source, target = edge.split('_')
        G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3, label=label)

    dot = nx.nx_pydot.to_pydot(G)
    dot.write_pdf(outfile)


# inputs
posterior_file = sys.argv[1]
primary_tissue = sys.argv[2]
n = int(sys.argv[3])

# # testing
# posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/metastabayes/MMUS1457/CP01/combined.trees"
# primary_tissue = "PRL"
# n = 10


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

# get n number of trees randomly from the posterior
mcc_top_n_trees = random.sample(trees, k=n) # a better approach would be to group the posterior by common features and then sample from those groups, but this is not yet implemented

# sample trees, convert to migration graphs, and output as files
i = 1
for tree in mcc_top_n_trees:
    outfile = posterior_file.split(".")[0] + "_migration_graph_" + str(i) + ".pdf"
    tree_to_migration_graph(tree, primary_tissue, outfile)
    i = i + 1
