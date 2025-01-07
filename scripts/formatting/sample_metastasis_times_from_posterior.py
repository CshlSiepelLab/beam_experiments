#!/usr/bin/env python3

import sys
import os
import re
from ete3 import Tree
from copy import deepcopy
import dendropy
from multiprocessing import Pool
from collections import deque
import csv
import pickle as pkl

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

    # remove quoted node names that is default in dendropy
    newick = newick.replace("\'", "")

    # the transition to ete3 tree here is just by preference
    # likely all the same subsequent analysis can be done with dendropy
    ete_tree = Tree(newick, format=3)

    return ete_tree

def level_order_traversal_met_events(tree, root_to_origin_height):

    met_times = {}
    migrations = set()

    for node in tree.traverse("levelorder"):
        if node.is_root():
            parent_tissue = origin_tissue
            parent_time = origin_time
            node_time = root_to_origin_height
        else:
            parent_node = node.up
            parent_tissue = parent_node.name.split("_")[-1]
            root = tree.get_tree_root()
            parent_time = root_to_origin_height + root.get_distance(parent_node.name)
            node_time = root_to_origin_height + root.get_distance(node.name)

        node_tissue = node.name.split("_")[-1]

        if node_tissue != parent_tissue:
            migration = f"{parent_tissue}_{node_tissue}"
            migration_time = (parent_time, node_time)
            if migration not in migrations:
                migrations.add(migration)
                migration = migration + "_1"
                met_times[migration] = migration_time
            else:
                existing_migrations = [key for key in met_times.keys() if migration in key]
                i = max([int(key.split("_")[-1]) for key in existing_migrations]) + 1
                migration = migration + "_" + str(i)
                met_times[migration] = migration_time
    
    return met_times



def is_ultrametric(tree):
    root = tree.get_tree_root()
    leaf_distances = set()
    for leaf in tree.iter_leaves():
        leaf_distances.add(round(root.get_distance(leaf.name), 3))
    return len(leaf_distances) == 1
        

# user inputs
posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/chain_3.trees"
origin_time = 250
origin_tissue = "P"
outfile = "./test_metastasis_times.pkl"
cores=1

# process beast posterior
burnin_percent = 0.1
beast_tree_list= dendropy.TreeList.get(path=posterior_file, schema="nexus")
num_beast_trees = len(beast_tree_list)
num_discard = round(num_beast_trees * burnin_percent)
beast_tree_list = beast_tree_list[num_discard:]
num_beast_trees = len(beast_tree_list)

all_met_events = {}

for tree in beast_tree_list:

    # put tissue location tree annotations in the node names and convert to ete3 tree (these are personal preferences)
    named_tree = rename_tree(tree)

    # verify that the tree is ultrametric, ie. same distance from the root to all leaves
    if not is_ultrametric(named_tree):
        raise ValueError("The tree is not ultrametric")

    # get the total tree height from root to any leaf
    tree_height = named_tree.get_farthest_leaf()[1]

    # calculate the time from the root to the origin, which is not output directly by BEAM so needs to be considered independently
    root_to_origin_height = origin_time - tree_height

    # traverse the tree to get the times of all metastatic envents in the migration graph implied by the tree
    timed_met_events = level_order_traversal_met_events(named_tree, root_to_origin_height)

    # record met events with the sample number in the posterior
    all_met_events[tree.label.split(" ")[1]] = timed_met_events

with open(outfile, "wb") as f:
    pkl.dump(all_met_events, f)




