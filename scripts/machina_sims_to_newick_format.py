#!/usr/bin/env python3

# This script takes in .tree and .vertex.labeling files from the simulted data in the MACHINA github repo and formats them to the newick files necessary to run BEAST and MACHINA in the existing pipeline setup

import sys
from ete3 import Tree
import pandas as pd

def build_tree_from_file(file_path):
    tree = Tree()
    edges = []
    letter_edges = []   # Solves edge case where MACHINA labeling is flipped?
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().replace('_', '-').split(' ')
            parent, child = fields[0], fields[1]
            if parent == "GL":
                gl_edge = (parent, child)
                continue
            if not parent[0].isdigit():
                letter_edge = (parent,child)
                letter_edges.append(letter_edge)
                continue
            edge = (parent, child)
            edges.append(edge)
    # Sort edges based on first value in key to prevent out of order additions to tree that result in duplicate nodes
    edges = sorted(edges, key=lambda x: int(x[0].split(";")[0]))
    edges.insert(0, gl_edge)
    letter_edges = sorted(letter_edges, key=lambda s: s[0][0])
    edges.extend(letter_edges)
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
        if parent.name != "GL" and len(parent.get_children()) == 1:
            leaf.detach()
            parent.name = leaf_name
    gl_node = tree.search_nodes(name="GL")[0]
    return gl_node, edges


tree_file = sys.argv[1]
vertex_labeling_file = sys.argv[2]

# tree_file = "machina_m8_sim_data/seed10046/T_seed10046.tree"
# vertex_labeling_file = "machina_m8_sim_data/seed10046/T_seed10046.vertex.labeling"

tissue_df = pd.read_csv(vertex_labeling_file, sep='\s+', names=['node', 'tissue'])
tree, edges = build_tree_from_file(tree_file)

for node in tree.traverse():
    node.name = node.name.replace(";","-")

unlabeled_names = [node.name for node in tree.traverse() if node.name != '']
tissue_df['node'] = tissue_df['node'].str.replace("_", "-")
tissue_df['node'] = tissue_df['node'].str.replace(";", "-")
tissue_df_subset = tissue_df[tissue_df['node'].isin(unlabeled_names)]
tissue_dict = dict(zip(tissue_df_subset['node'], tissue_df_subset['tissue']))

labeled_tree = tree.copy()
for node in labeled_tree.traverse():
    name = node.name
    node.name = name + "_" + tissue_dict[name]


output_prefix = tree_file.split(".")[0]
# Output unlabeled tree to newick
output_unlabeled_file = output_prefix + "_unlabeled_true_tree.nwk"
tree.write(outfile=output_unlabeled_file, format=5)

# Output tissue label tsv
output_tissues_tsv = output_prefix + "_tissues.tsv"
tissue_df_subset.to_csv(output_tissues_tsv, sep='\t', index=False)
    
# Output fully tissue labeled newick
output_labeled_file = output_prefix + "_tissue_labeled_true_tree.nwk"
labeled_tree.write(outfile=output_labeled_file, format=8)


