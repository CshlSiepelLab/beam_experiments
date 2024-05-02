#!/usr/bin/env python3

### This script takes in a newick tree and tsv style tissue map for the tips of the tree and then outputs inferred newicks for both a consensus and random approach to labeling tissues for internal nodes

import sys
import ete3
import pandas as pd
from collections import Counter
import random

def label_tissues_consensus(tree, tissues_df):
    copy_tree = tree.copy()
    for node in copy_tree.traverse():
        # leaves are assumed to be known labels
        if node.is_leaf():
            tissue = tissues_df.loc[tissues_df['cell'] == str(node.name), 'tissue'].values[0]
            node.name = f"{node.name}_{tissue}"
        # internal nodes are chosen by consensus of leaf tissues
        else:
            node_name = node.name
            children = [tissues_df.loc[tissues_df['cell'] == str(leaf.name), 'tissue'].values[0] for leaf in node.get_leaves()]
            counts = Counter(children)
            most_common_elements = counts.most_common(2)
            # if there is a clear winner, choose that tissue
            if len(most_common_elements) == 1 or most_common_elements[0][1] != most_common_elements[1][1]:
                consensus_tissue = most_common_elements[0][0]
            else:
                # tie breaker goes to primary
                if any('P' in elem for elem in most_common_elements):
                    consensus_tissue = 'P'
                # if tie doesn't involve primary, choose randomly
                else:
                    tied_elements = [elem[0] for elem in most_common_elements]
                    consensus_tissue = random.choice(tied_elements)
            # rename with consensus tissue label
            node.name = f"{node_name}_{consensus_tissue}"
    return copy_tree

def label_tissues_random(tree, tissues_df):
    copy_tree = tree.copy()
    tissues = tissues_df['tissue'].unique().tolist()
    for node in copy_tree.traverse():
        # leaves are assumed to be known labels
        if node.is_leaf():
            tissue = tissues_df.loc[tissues_df['cell'] == str(node.name), 'tissue'].values[0]
            node.name = f"{node.name}_{tissue}"
        # internal nodes are chosen randomly
        else:
            random_tissue = random.choice(tissues)
            node.name = f"{node.name}_{random_tissue}"
    return copy_tree

# tree_file = sys.argv[1]
# leaf_tissues_tsv = sys.argv[1]

tree_file = "cell_tree_seed7518.nwk"
leaf_tissues_tsv = "cell_tree_seed7518.labeling"


tree = ete3.Tree(tree_file, format=3)
tissue_map = pd.read_csv(leaf_tissues_tsv, sep='\s+', header=None, names=['cell', 'tissue'], dtype={'cell': str, 'tissue': str})

# get results
random_tree = label_tissues_random(tree, tissue_map)
consensus_tree = label_tissues_consensus(tree, tissue_map)

# output trees to newick files
random_output = tree_file.split(".nwk")[0] + "_random_tissues.nwk"
random_tree.write(outfile = random_output, format=8, format_root_node=True)

consensus_output = tree_file.split(".nwk")[0] + "_consensus_tissues.nwk"
consensus_tree.write(outfile = consensus_output, format=8, format_root_node=True)
        