#!/usr/bin/env python3

import sys
from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, TextFace
import pandas as pd
import dendropy
import os


DEFAULT_COLORS = ["#006400", "#FF0000", "#0000CD", "#FFA500", "#800080", "#808080", "#FFC0CB", "#ADD8E6", "#A52A2A", "#FFFF00"]*3


def plot_tree_and_graph(newick_file, tissues_file, primary_tissue, total_time, outfile):

    tree = dendropy.Tree.get(path=newick_file, schema="newick")

    # Replace any '.' or '-' characters in node names with '_'
    for node in tree.preorder_node_iter():
        if node.label and ('.' in node.label or '-' in node.label):
            node.label = node.label.replace('.', '_').replace('-', '_')
    
    # Convert dendropy tree to ete3 tree
    newick_str = tree.as_string(schema="newick")
    tree = Tree(newick_str, format=1)

    tissues_dict = {}
    with open(tissues_file, 'r') as f:
        for line in f:
            name, tissue = line.strip().split(",")
            tissues_dict[name] = tissue

    all_tissues = set()
    all_tissues.add(primary_tissue)

    # Get tissue labels and assign node names
    i = 1
    for node in tree.traverse():
        node.tissue = None  # Initialize tissue attribute
        if node.name in tissues_dict:
            node.tissue = tissues_dict[node.name]
            if node.tissue not in all_tissues:
                all_tissues.add(node.tissue)
        
        if node.name == "":
            if node.is_root():
                node.name = "root"
                if "0" in tissues_dict:
                    node.tissue = tissues_dict["0"]
                else:
                    node.tissue = primary_tissue
            else:
                node.name = f"node{i}"
                i += 1
    
    # Check that the tree is ultrametric
    dists = set()
    for node in tree.traverse():
        if node.is_leaf():
            dists.add(round(node.get_distance(tree), 3))
    if len(dists) != 1:
        print("WARNING: Tree sample is not ultrametric.")

    # Add origin node above the root
    tree_height = dists.pop()
    origin = Tree(name="origin", dist=0)
    origin.tissue = primary_tissue
    root = tree.get_tree_root()
    root.dist = total_time - tree_height
    origin.add_child(root)
        
    # get all tissue names and assign them colors
    all_tissues = sorted(list(set(all_tissues) - {primary_tissue}))
    custom_colors = {node: color for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)]) if node != primary_tissue}
    all_tissues = [primary_tissue] + all_tissues
    custom_colors[primary_tissue] = "black"

    # Plot tree
    ts = TreeStyle()
    ts.rotation = 90
    ts.scale = 1
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_border = False
    ts.show_scale = False
    ts.mode = "r"

    # Add legend for each color
    for tissue in all_tissues:
        ts.legend.add_face(CircleFace(10, custom_colors[tissue]), column=0)
        ts.legend.add_face(TextFace(tissue, fsize=12), column=1)

    # Setup node style
    all_labeled = True
    for node in origin.traverse():
        nstyle = NodeStyle()
        nstyle["shape"] = "circle"
        nstyle["size"] = 10
        if node.tissue == None:
            color = "grey"
            all_labeled = False
        else:
            color = custom_colors[node.tissue]
        nstyle["fgcolor"] = color
        # if not all nodes are labeled than color edges grey since the tree must not have been fully labeled yet
        if all_labeled == False:
            color = "grey"
        nstyle["hz_line_color"] = color
        nstyle["vt_line_color"] = color
        nstyle["hz_line_width"] = 3
        nstyle["vt_line_width"] = 3
        
        node.set_style(nstyle)
        
    # Set the QT_QPA_PLATFORM environment variable
    os.environ["QT_QPA_PLATFORM"] = "offscreen"

    origin.render(outfile, tree_style=ts)


# inputs
# newick_file = sys.argv[1]
# tissues_file = sys.argv[2]
# primary_tissue = sys.argv[3]
# total_time = int(sys.argv[4])
# outfile = sys.argv[5]

newick_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/mach2/5k/58/M-T-0.nwk"
tissues_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/mach2/5k/58/M-T-0_labeling.csv"
primary_tissue = "LL"
total_time = 54
outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/mach2/5k/58/M-T-0_labeling.pdf"


plot_tree_and_graph(newick_file, tissues_file, primary_tissue, total_time, outfile)
