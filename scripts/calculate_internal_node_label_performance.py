#!/usr/bin/env python3

# This script takes in trees for the prediction results from BEAST2 FixedTreeAnalysis maximum clade credibility tree from TreeAnotator and compares accuracy to MACHINA predictions against the ground truth simulated tree.

import sys
from ete3 import Tree
import dendropy

def get_labels_newick(tree):
    labels = []
    for node in tree.traverse():
        if node.is_leaf() or node.is_root():
            continue
        else:
            labels.append(node.name)
    return labels

def get_posterior_labels_strict(tree):
    # Takes in nexus tree with posterior probabilities for tissue locations and returns list of most probable locations concatenated to node name
    labels = []
    for node in tree.internal_nodes():
            prediction = node.label + "_" + node.annotations.get_value('location')
            labels.append(prediction)
    return labels

def get_posterior_labels_relaxed(tree,threshold):
    # Takes in nexus tree with posterior probabilities for tissue locations and returns list of most probable locations while allowing more than one most probable based on a posterior probability difference threshold where results are then concatenated to node name
    # Threshold is a floating point indicating the range of posterior probabilities from the max for which to accept more than one correct prediction
    labels = []
    for node in tree.internal_nodes():
        location_prob = node.annotations.get_value('location.prob')
        if float(location_prob) != 1.0:
            location_set = node.annotations.get_value('location.set')
            location_probs = node.annotations.get_value('location.set.prob')
            location_probs = [float(prob) for prob in location_probs]
            max_prob = max(location_probs)
            min_prob_threshold = max_prob - threshold
            max_tissues = [location_set[i] for i, prob in enumerate(location_probs) if prob > min_prob_threshold]
            current_label = str(node.label)
            for tissue in max_tissues:
                prediction = current_label + f"_{tissue}"
                labels.append(prediction)
        else:
            tissue = node.annotations.get_value('location')
            current_label = str(node.label)
            prediction = current_label + f"_{tissue}"
            labels.append(prediction)
    return labels


true_file=sys.argv[1]
beast_file=sys.argv[2]
machina_file=sys.argv[3]

# true_file="compare_beast_machina_fixedtree_2_2_24/sim_results_sim1/sim1_true_tissues.nwk"
# beast_file="compare_beast_machina_fixedtree_2_2_24/sim_results_sim1/tissue_tree_with_trait.tree"
# machina_file="compare_beast_machina_fixedtree_2_2_24/sim_results_sim1/machina_tree_all_tissue_labels.nwk"

data_id = true_file.split("/")[-1].split(".")[0]

true_tree = Tree(true_file, format=8)
beast_tree = dendropy.Tree.get(path=beast_file, schema='nexus')
machina_tree = Tree(machina_file, format=8)

# Relabel beast tree internal nodes based on true tree dictionary
true_tree_no_tissues = true_tree.copy()
node_leaf_dict = {}
for node in true_tree_no_tissues.traverse():
    name = node.name.split("_")[0]
    if node.is_leaf() == False:
        leaves = "_".join(sorted([leaf.name.split("_")[0] for leaf in node.get_leaves()]))
        node_leaf_dict[leaves] = name

for node in beast_tree.internal_nodes():
    leaves = "_".join(sorted([leaf.taxon.label.split("ll")[1] for leaf in node.leaf_nodes()]))
    node.label = node_leaf_dict[leaves]

# Get internal node tissue labels for true tree and machina
true_labels = get_labels_newick(true_tree)
machina_labels = get_labels_newick(machina_tree)

total = len(true_labels)

machina_correct = [label for label in true_labels if label in machina_labels]
machina_accuracy = len(machina_correct) / total

beast_strict_labels = get_posterior_labels_strict(beast_tree)
beast_relaxed_labels = get_posterior_labels_relaxed(beast_tree, 0.05)

beast_strict_correct = [label for label in true_labels if label in beast_strict_labels]
beast_relaxed_correct = [label for label in true_labels if label in beast_relaxed_labels]

beast_strict_accuracy = len(beast_strict_correct) / total
beast_relaxed_accuracy = len(beast_relaxed_correct) / total

# Compute non-primary (np) tissue (t1) accuracy
np_true_labels = [label for label in true_labels if 't1' not in label]
np_total = len(np_true_labels)
np_machina_labels = [label for label in machina_labels if label in np_true_labels]
np_beast_strict_labels = [label for label in beast_strict_labels if label in np_true_labels]
np_beast_relaxed_labels = [label for label in beast_relaxed_labels if label in np_true_labels]

np_machina_accuracy = len(np_machina_labels) / np_total
np_beast_strict_accuracy = len(np_beast_strict_labels) / np_total
np_beast_relaxed_accuracy = len(np_beast_relaxed_labels) / np_total

outputfile = "/".join(true_file.split("/")[:-1]) + "/compare_machina_beast_internal_node_performance.tsv"
with open(outputfile, "w") as file:
    header_str = "data_id\tmachina\tbeast_strict\tbeast_relaxed\tmachina_nonprimary\tbeast_strict_nonprimary\tbeast_relaxed_nonprimary"
    accuracy_str = f"{data_id}\t{machina_accuracy}\t{beast_strict_accuracy}\t{beast_relaxed_accuracy}\t{np_machina_accuracy}\t{np_beast_strict_accuracy}\t{np_beast_relaxed_accuracy}"
    file.write(f"{header_str}\n{accuracy_str}")