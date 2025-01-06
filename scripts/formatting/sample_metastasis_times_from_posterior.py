#!/usr/bin/env python3

import sys
import os
import re
from ete3 import Tree
from copy import deepcopy
import dendropy
from multiprocessing import Pool
from collections import deque

def rename_tree(tree):
    tree_copy = deepcopy(tree)
    i = 0
    for node in tree_copy.preorder_node_iter():
        try:
            prediction = node.taxon.label + "_" + node.annotations.get_value('location')
            node.taxon.label = prediction
        except Exception as e:
            prediction = f"node{i}" + "_" + node.annotations.get_value('location')
            i += 1
        node.label = prediction
    
    newick = tree_copy.as_string(schema='newick', suppress_edge_lengths=False, node_label_element_separator=',')

    # remove quoted node names
    newick = newick.replace("\'", "")

    ete_tree = Tree(newick, format=3)

    return ete_tree

def level_order_traversal(tree):
    events = {}
    queue = deque([tree.seed_node])
    while queue:
        node = queue.popleft()
        for child in node.child_node_iter():
            edge = f"{node.label.split('_')[-1]}_{child.label.split('_')[-1]}"
            if edge not in events:
                events[edge] = [child.edge_length]
            else:
                events[edge].append(child.edge_length)
            queue.append(child)

def is_ultrametric(tree):
    root = tree.seed_node
    leaf_distances = set()
    for leaf in tree.leaf_node_iter():
        leaf_distances.add(tree.distance_from_root(leaf))
    return len(leaf_distances) == 1

def tree_height(tree):
    root = tree.seed_node
    max_distance = 0
    for leaf in tree.leaf_node_iter():
        distance = tree.distance_from_root(leaf)
        if distance > max_distance:
            max_distance = distance
    return max_distance
        

# user inputs
posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_12_31_24_uniform_50cells_50sites_data_7_24_24/beam/mS_854/chain_3.trees"
origin_time = 250
primary_tissue = "P"
outfile = "./test_metastasis_times.csv"
cores=1

# process beast posterior
burnin_percent = 0.1
beast_tree_list= dendropy.TreeList.get(path=posterior_file, schema="nexus")
num_beast_trees = len(beast_tree_list)
num_discard = round(num_beast_trees * burnin_percent)
beast_tree_list = beast_tree_list[num_discard:]
num_beast_trees = len(beast_tree_list)

tree = beast_tree_list[0]

named_tree = rename_tree(tree)

if not is_ultrametric(tree):
    raise ValueError("The tree is not ultrametric")

total_height = tree_height(tree)
print(f"Total height from root to leaves: {total_height}")

named_tree = rename_tree(tree)

timed_met_events = level_order_traversal(named_tree)




