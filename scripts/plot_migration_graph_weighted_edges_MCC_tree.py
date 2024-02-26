#!/usr/bin/env python3

### This script takes in a MCC tree from TreeAnnotator made from the BEAST2 posterior and then collapses the tree to a migration graph with edges weighted based on node tissue location probabilities. We essentially make a graph of all possible routes and then each route intensity is based on probability of occuring in the tree given node probabilities.

import re, sys
import pandas as pd

def remove_bracket_content(match):
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

# consensus_tree_file = sys.argv[1]

consensus_tree_file = "beast_gundem_2015_2_21_24/A10_sym/tissue_tree_with_trait.tree"

with open(consensus_tree_file, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('tree'):
            tree_info = line

# remove tree name and = from tree info to get only the newick
tree = ''.join(tree_info.split(' ')[3:])

# strip tree string to newick with associated dataframe of annotations
bracket_content_pattern = re.compile(r'\[.*?\]')
annotations = []
newick = re.sub(bracket_content_pattern, remove_bracket_content, tree)
annotations = [re.split(r',(?![^{]*})', x.replace("&", "")) for x in annotations]

annotations = [{key: value.replace("{", "").replace("}", "") for trait in annotation for key, value in [trait.split("=")]} for annotation in annotations]

# label nodes in newick with only leaf labels
node_labeled_newick = label_nodes(newick)
node_labels = get_labels(node_labeled_newick)

# make a dictionary for annotations to node labels
annotations_dict = {}
for node in node_labels:
    annotations_dict[node] = annotations[node_labels.index(node)]

### NEED TO PLOT GRAPH HERE NOW THAT ANNOTATIONS ARE SAVED IN A DICT AND NEWICK IS ITERABLE
