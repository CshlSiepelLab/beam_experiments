#!/usr/bin/env python3

# This script takes in .tree and .vertex.labeling files from the real data in the MACHINA github repo and formats them to a leaf labeled newick file

import sys
from ete3 import Tree
import pandas as pd

def build_tree_from_file(file_path):
    edges = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().replace('_', '-').split(' ')
            parent, child = fields[0], fields[1]
            edge = (parent, child)
            edges.append(edge)
    # assume first parent is the root
    root = edges[0][0]
    tree = Tree(name=root)
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
    # Remove repetitive leaf labels for leaves coming form a parent with only one child
    for leaf in tree.iter_leaves():
        leaf_name = leaf.name
        parent = leaf.up
        if parent.name != root and len(parent.get_children()) == 1:
            leaf.detach()
            parent.name = leaf_name
    # Remove polytomies by forcing bifurcations with branch length of 0
    for node in tree.traverse("postorder"):
        if len(node.children) > 2:
            # Create virtual nodes with zero-length branches
            virtual_node = Tree()
            while len(node.children) > 1:
                child = node.children.pop()
                virtual_node.add_child(child)
            node.add_child(virtual_node, dist=0)
    return tree, edges


tree_file = sys.argv[1]
vertex_labeling_file = sys.argv[2]
primary_tissue = sys.argv[3]

# tree_file = "gundem_a10/A10.tree"
# vertex_labeling_file = "gundem_a10/A10.labeling"
# primary_tissue = "prostate"

tissue_df = pd.read_csv(vertex_labeling_file, sep='\s+', names=['node', 'tissue'])
tree, edges = build_tree_from_file(tree_file)

tissue_df.loc[: ,'node'] = [s.replace("_", "-") for s in list(tissue_df.loc[: ,'node'].values)]
tissue_df.loc[: ,'tissue'] = [s.replace("_", "-") for s in list(tissue_df.loc[: ,'tissue'].values)]

tissue_dict = dict(zip(tissue_df['node'], tissue_df['tissue']))

labeled_tree = tree.copy()
for leaf in labeled_tree.iter_leaves():
    name = leaf.name
    tissue = tissue_dict[name]
    leaf.name = name + "_" + tissue


output_prefix = tree_file.split(".")[0]
# Output unlabeled tree to newick
output_unlabeled_file = output_prefix + "_unlabeled_tree.nwk"
tree.write(outfile=output_unlabeled_file, format=5)

# Output fully tissue labeled newick
output_labeled_file = output_prefix + "_tissue_labeled_tree.nwk"
labeled_tree.write(outfile=output_labeled_file, format=5)

# Output tissue label tsv
output_tissues_tsv = output_prefix + "_tissues.tsv"
tissue_df.to_csv(output_tissues_tsv, sep='\t', index=False)


