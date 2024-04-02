#!/usr/bin/env python3

import sys
from ete3 import Tree
import dendropy
from copy import deepcopy

def remove_zero_length_nodes(tree):
    for node in tree.internal_nodes():
        if node.edge_length == 0:
            parent = node.parent_node
            if parent is not None:
                parent.remove_child(node)
                children = node.child_nodes()
                for child in children:
                    parent.add_child(child)

def dendropy_beast_to_ete_newick_with_strict_locations(tree):
    tree_copy = deepcopy(tree)
    i = 0
    for node in tree_copy.preorder_node_iter():
        try:
            prediction = node.taxon.label + "_" + node.annotations.get_value('location')
        except Exception as e:
            prediction = f"node{i}" + "_" + node.annotations.get_value('location')
            i += 1
        node.label = prediction
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("\'", ""), format=3)
    return ete_tree

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

def calculate_metrics(true_counts, inferred_counts):
        # get general totals for formulas
    total_inferred = sum(inferred_counts.values())
    total_true = sum(true_counts.values())

    # get union of positives for true and inferred
    union_positives = 0
    for key in inferred_counts.keys():
        inferred_value = inferred_counts[key]
        if key not in true_counts:
            continue
        else:
            true_value = true_counts[key]
            if true_value < inferred_value:
                union_positives += true_value
            else:
                union_positives += inferred_value
    # calculate recall (union of inferred and true positives over inferred total positives)
    recall = union_positives / total_inferred
    # calculate precision (union of inferred and true positives over total true positives)
    precision = union_positives / total_true
    # calculate F1 score (2((precision * recall)/(precision + recall)))
    f1 = 2 * ((precision * recall) / (precision + recall))
    return f1, recall, precision

def main():
    true_tree_file=sys.argv[1]
    beast_trees_file=sys.argv[2]

    # true_tree_file="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874_tissue_labeled_tree.nwk"
    # beast_trees_file="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/joint_inference_beast_tissues.trees"

    burnin_percent=0.1
    primary_tissue="P"

    # true counts are the same for all beast tree comparisons
    true_tree = Tree(true_tree_file, format=8)
    true_tree.get_tree_root().name = f'0_{primary_tissue}'
    true_counts=get_migration_counts(true_tree)

    beast_tree_list = dendropy.TreeList()
    beast_tree_list.read(path=beast_trees_file, schema="nexus")

    num_beast_trees = len(beast_tree_list)
    num_discard = round(num_beast_trees * burnin_percent)
    beast_tree_list = beast_tree_list[num_discard:]

    posteriors = []
    f1_scores=[]

    for tree in beast_tree_list:
        posteriors.append(float(tree.annotations.get_value('posterior')))
        remove_zero_length_nodes(tree)
        inferred_tree = dendropy_beast_to_ete_newick_with_strict_locations(tree)

        # calculate beast counts for each tree, get metrics, and weight by posterior probability
        inferred_tree.get_tree_root().name = f'0_{primary_tissue}'
        inferred_counts=get_migration_counts(inferred_tree)

        f1, recall, precision = calculate_metrics(true_counts, inferred_counts)
        f1_scores.append(f1)

    avg_posterior_f1 = sum(f1_scores) / len(f1_scores)
    print(f"F1 score: {avg_posterior_f1}")

if __name__ == "__main__":
    main()
