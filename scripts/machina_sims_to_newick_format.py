#!/usr/bin/env python3

# This script takes in .tree and .vertex.labeling files from the simulted data in the MACHINA github repo and formats them to the newick files necessary to run BEAST and MACHINA in the existing pipeline setup

from ete3 import Tree
import pandas as pd

def build_tree_from_file(file_path):
    tree = Tree()
    edges = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split(' ')
            parent, child = fields[0], fields[1]
            if parent == "GL":
                continue
            edge = (parent, child)
            edges.append(edge)
    # Sort edges based on first value in key to prevent out of order additions to tree that result in duplicate nodes
    edges = sorted(edges, key=lambda x: int(x[0].split(";")[0]))
    for edge in edges:
        parent, child = edge
        parent = str(parent)
        child = str(child)
        parent_node = tree.search_nodes(name=parent)
        child_node = tree.search_nodes(name=child)
        if not parent_node:
            parent_node = tree.add_child(name=parent)
        else:
            parent_node = parent_node[0]
        if not child_node:
            child_node = parent_node.add_child(name=child)
    # Remove repetitive leaf labels for leaves named such as "13_P" whose parent is "13"
    for leaf in tree.iter_leaves():
        leaf_name = leaf.name
        parent = leaf.up
        if len(parent.get_children()) == 1:
            leaf.detach()
            parent.name = leaf_name
    return tree, edges


tree_file = "machina_m5_sim_data/T_seed0.tree"
vertex_labeling_file = "machina_m5_sim_data/T_seed0.vertex.labeling"

tissue_labels = pd.read_csv(vertex_labeling_file, sep='\s+', names=['node', 'tissue'])
tree, edges = build_tree_from_file(tree_file)


tree_nodes = [node.name for node in tree.traverse()]
label_names = list(tissue_labels['node'].values)

