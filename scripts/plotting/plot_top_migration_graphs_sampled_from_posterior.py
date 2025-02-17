#!/usr/bin/env python3

import re
import sys
import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, TextFace
import pandas as pd
import matplotlib
import random
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import to_pydot


DEFAULT_COLORS = ["#006400", "#FF0000", "#0000CD", "#FFA500", "#800080", "#808080", "#FFC0CB", "#ADD8E6", "#A52A2A", "#FFFF00"]*3

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

def plot_tree_and_graph(tree, primary_tissue, total_time, outprefix, num):
    newick_str = tree.split(' ')[3]
    tree = Tree(newick_str, format=1)
    all_tissues = set()
    all_tissues.add(primary_tissue)

    # Get tissue labels and assign node names
    i = 1
    for node in tree.traverse():
        node.name, annotation = node.name.split("[")
        node.tissue = re.search(r'&location="([^"]+)"', annotation).group(1)
        
        if node.name == "":
            if node.is_root():
                node.name = "root"
            else:
                node.name = f"node{i}"
                i += 1
        # # to use the original names (commented out here since it can cause plotting issues)
        # else:
        #     node.name = names_dict[node.name]
                
        if node.tissue not in all_tissues:
            all_tissues.add(node.tissue)
    
    # Check that the tree is ultrametric
    dists = set()
    for node in tree.traverse():
        if node.is_leaf():
            dists.add(round(node.get_distance(tree), 5))
    if len(dists) != 1:
        print("WARNING: Tree sample is not ultrametric when it should be for BEAST2 output. Check the newick string to verify: ", newick_str)

    # Add origin node above the root
    tree_height = dists.pop()
    origin = Tree(name="origin", dist=0)
    origin.tissue = primary_tissue
    root = tree.get_tree_root()
    root.dist = total_time - tree_height
    origin.add_child(root)

    # traverse the tree for the migration graph
    migration_counts = {}

    for node in origin.traverse():
        # root is now the origin so skip it here
        if node.is_root():
            continue
        if node.up.tissue != node.tissue:
            migration = f"{node.up.tissue}_{node.tissue}"
            if migration in migration_counts:
                migration_counts[migration] += 1
            else:
                migration_counts[migration] = 1
        
    # get all tissue names and assign them colors
    all_tissues = sorted(list(set(all_tissues) - {primary_tissue}))
    custom_colors = {node: color for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)]) if node != primary_tissue}
    all_tissues = [primary_tissue] + all_tissues
    custom_colors[primary_tissue] = "black"

    # Plot migration graph
    G = nx.MultiDiGraph()

    for node in all_tissues:
        G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0, fontsize=32)

    for edge, count in migration_counts.items():
        if count == 1:
            label = ""
        else:
            label = str(count)
        source, target = edge.split('_')
        G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3, label=label, fontsize=24)

    dot = nx.nx_pydot.to_pydot(G)
    outfile = outprefix + "_migration_graph_" + str(num) + ".pdf"
    dot.write_pdf(outfile)

    # Plot tree
    ts = TreeStyle()
    ts.rotation = 90
    ts.scale = 1
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_border = False
    ts.show_scale = False
    ts.mode = "r"

    # Add legend for each color
    for tissue in all_tissues:
        ts.legend.add_face(CircleFace(10, custom_colors[tissue]), column=0)
        ts.legend.add_face(TextFace(tissue, fsize=12), column=1)

    # Setup node style
    for node in origin.traverse():
        nstyle = NodeStyle()
        nstyle["shape"] = "circle"
        nstyle["size"] = 10
        nstyle["hz_line_color"] = custom_colors[node.tissue]
        nstyle["vt_line_color"] = custom_colors[node.tissue]
        nstyle["hz_line_width"] = 3
        nstyle["vt_line_width"] = 3
        nstyle["fgcolor"] = custom_colors[node.tissue]
        node.set_style(nstyle)
        
    outfile = outprefix + "_tree_" + str(num) + ".pdf"
    origin.render(outfile, tree_style=ts)

    # Plot metastasis timing
    metastasis_times = {}

    for node in origin.traverse():
        # root is now the origin so skip it here
        if node.is_root():
            continue
        if node.up.tissue != node.tissue:
            migration = f"{node.up.tissue}_{node.tissue}"
            # time of metastasis is halfway on the branch to the node in a new site
            time = node.up.get_distance(origin) + (node.dist / 2)
            if migration in metastasis_times:
                metastasis_times[migration].append(time)
            else:
                metastasis_times[migration] = [time]
    
    # Plot rectangle spike plot for migration counts
    fig, ax = plt.subplots(figsize=(12, 2))

    for migration, times in metastasis_times.items():
        source, target = migration.split('_')
        color = custom_colors[target]
        for time in times:
            ax.plot([time, time], [0, 1], color=color, linewidth=4)

    fs = 18
    ax.set_xlim(0, total_time)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time', fontsize=fs)
    ax.set_ylabel('Migrations', fontsize=fs)
    ax.yaxis.set_visible(False)
    ax.tick_params(axis='x', labelsize=fs)

    # Add legend
    handles = [plt.Line2D([0], [0], color=color, lw=4) for color in custom_colors.values()]
    labels = custom_colors.keys()
    ax.legend(handles, labels, title='', title_fontsize=fs, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., frameon=False, fontsize=fs)

    plt.tight_layout()
    outfile = outprefix + "_migration_timing_" + str(num) + ".pdf"
    plt.savefig(outfile)
    plt.close()



# inputs
posterior_file = sys.argv[1]
primary_tissue = sys.argv[2]
total_time = int(sys.argv[3])
n = int(sys.argv[4])

# # testing
# posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_1_22_25/beam_gtr/5k/58/combined.trees"
# primary_tissue = "LL"
# total_time = 54
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

# # discard 10% of tree for burnin
# burnin = int(len(trees) * 0.1)
# trees = trees[burnin:]

# get n number of trees randomly from the posterior
mcc_top_n_trees = random.sample(trees, k=n) # a better approach would be to group the posterior by common features and then sample from those groups, but this is not yet implemented

# sample trees, convert to migration graphs, and output as files
i = 1
for tree in mcc_top_n_trees:
    plot_tree_and_graph(tree, primary_tissue, total_time, posterior_file.split(".")[0], i)
    i = i + 1
