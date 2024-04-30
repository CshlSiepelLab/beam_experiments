#!/usr/bin/env python3

### This script takes in a newick tree and tsv style tissue map for the tips of the tree and then outputs inferred newicks for both a consensus and random approach to labeling tissues for internal nodes

import sys
import ete3
import pandas as pd
from collections import Counter
import random

def label_tissues_consensus(tree):


def label_tissues_random(tree):


# tree_file = sys.argv[1]
# leaf_tissues_tsv = sys.argv[1]

tree_file = "cell_tree_seed7518.nwk"
leaf_tissues_tsv = "cell_tree_seed7518.labeling"


tree = ete3.Tree(tree_file, format=3)
tissue_map = pd.read_csv(leaf_tissues_tsv, sep="\t", header=None)



for node in tree.traverse():
    if node.is_leaf():
        pass
    else:
        node_name = node.name
        children = [leaf.name.split("_")[1] for leaf in node.get_leaves()]
        counts = Counter(children)
        most_common_elements = counts.most_common(2)
        if len(most_common_elements) > 1 and most_common_elements[0][1] == most_common_elements[1][1]:
            t1_in_common_elements = any('t1' in elem for elem in most_common_elements)
            if t1_in_common_elements:
                consensus_tissue = 't1'
            else:
                tied_elements = [elem[0] for elem in most_common_elements]
                consensus_tissue = random.choice(tied_elements)
        else:
            consensus_tissue = most_common_elements[0][0]
        new_node_name = f"{node_name}_{consensus_tissue}"
        node.name = new_node_name



tree.write(format=8, outfile = f'{output_dir}/consensus_tissue_tree_all_tissue_labels.nwk')
        