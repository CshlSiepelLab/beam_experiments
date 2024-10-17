#!/usr/bin/env python3

import sys
import os
import re
from scipy.stats import gaussian_kde
from ete3 import Tree
from copy import deepcopy
import dendropy
from multiprocessing import Pool

def get_migration_counts(tree):
    migration_counts = {}
    for node in tree.traverse():
        if node.is_root():
            continue
        else:
            node_tissue = node.name.split("_")[-1]
            parent_tissue = node.up.name.split("_")[-1]
            if node_tissue == parent_tissue:
                continue
            migration = f"{parent_tissue}_{node_tissue}"
            if migration not in migration_counts:
                migration_counts[migration] = 1
            else:
                migration_counts[migration] += 1
    return migration_counts

def optionally_add_origin_migration(tree, counts, primary):
    """
    Adds the origin branch to the migration counts if the root node location is different from the known origin primary tissue.

    Parameters:
    tree (ete3.Tree): The phylogenetic tree object.
    counts (dict): A dictionary containing migration counts.
    primary (str): The known primary tissue for the origin.

    Returns:
    dict: Updated migration counts including the origin branch if applicable.
    """
    counts_copy = counts
    root_node_location = tree.get_tree_root().name.split("_")[-1]
    if root_node_location != primary:
        migration = f"{primary}_{root_node_location}"
        if migration not in counts_copy:
            counts_copy[migration] = 1
        else:
            counts_copy[migration] += 1
    return counts_copy

def remove_zero_length_nodes(tree):
    for node in tree.internal_nodes():
        if node.edge_length == 0:
            parent = node.parent_node
            if parent is not None:  # prevents root from being lost
                parent.remove_child(node)
                children = node.child_nodes()
                for child in children:
                    parent.add_child(child)

def dendropy_beast_to_ete_newick_with_strict_locations(tree):
    # Note: this method does not convert the branch lengths properly since they are not used
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
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("\'", ""), format=3)
    return ete_tree


def get_consensus_graph(all_inferred_counts):
    # calculate total counts weighted by posterior probability which is uniform across all samples due to the mcmc sampling
    prob = 1 / len(all_inferred_counts)
    total_counts = {}
    for inferred_counts in all_inferred_counts:
        for pattern, count in inferred_counts.items():
            for num in range(1, count+1):
                edge = f"{pattern}_{num}"
                if edge not in total_counts:
                    total_counts[edge] = prob
                else:
                    total_counts[edge] += prob

    return total_counts

def process_tree_parallel(args):
    tree, primary_tissue = args
    posterior = float(tree.annotations.get_value('posterior'))
    remove_zero_length_nodes(tree)
    inferred_tree = dendropy_beast_to_ete_newick_with_strict_locations(tree)
    inferred_counts = get_migration_counts(inferred_tree)
    inferred_counts = optionally_add_origin_migration(inferred_tree, inferred_counts, primary_tissue)
    return posterior, inferred_counts

# user inputs
beast_posterior_file = sys.argv[1]
primary_tissue=sys.argv[2]
outfile = sys.argv[3]
cores = int(sys.argv[4])

# process beast posterior
burnin_percent = 0.1
beast_tree_list= dendropy.TreeList.get(path=beast_posterior_file, schema="nexus")
num_beast_trees = len(beast_tree_list)
num_discard = round(num_beast_trees * burnin_percent)
beast_tree_list = beast_tree_list[num_discard:]
num_beast_trees = len(beast_tree_list)

# create a pool of worker processes to process the trees
pool = Pool(processes=cores)
output = pool.map(process_tree_parallel, [(tree, primary_tissue) for tree in beast_tree_list])
pool.close()
pool.join()

posteriors, all_inferred_counts = zip(*output)


# fit a gaussian kernel density estimate to the posterior values to get a probability density function
posterior_prob_graph = get_consensus_graph(all_inferred_counts)

# output posterior_prob_graph to a file
with open(outfile, "w") as file:
    posterior_prob_graph = dict(sorted(posterior_prob_graph.items(), key=lambda x: x[1], reverse=True))
    for key, value in posterior_prob_graph.items():
        file.write(f"{key},{value}\n")
