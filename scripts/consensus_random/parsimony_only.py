#!/usr/bin/env python3

import sys
import ete3
import pandas as pd
import random

def label_tissues_parsimony(tree, tissues_df):
    '''
    Fitch parsimony algorithm to infer ancestral states of internal nodes for a tree
    '''
    def postorder(node):
        if node.is_leaf():
            # Assign known tissue type from the tissues_df to the leaf node
            node.final_tissue = tissues_df.loc[tissues_df['cell'] == str(node.name), 'tissue'].values[0]
            node.tissue_set = {node.final_tissue}
            node.name = f"{node.name}_{node.final_tissue}"
        else:
            # Process all children
            children_tissue_sets = [postorder(child) for child in node.children]
            
            # Compute the possible tissues for internal nodes based on the children's tissue sets
            intersection = set.intersection(*children_tissue_sets)
            if intersection:
                node.tissue_set = intersection
            else:
                node.tissue_set = set.union(*children_tissue_sets)
        return node.tissue_set
    
    def preorder(node, parent_tissue=None):
        # Leaf tissues are already known so skip them for tissue assignment
        if not node.is_leaf():
            if node.is_root():
                # The root tissue is known
                node.final_tissue = f"{primary_tissue}"
            elif parent_tissue and parent_tissue in node.tissue_set:
                # If parent tissue is in the node's set, choose it
                node.final_tissue = parent_tissue
            else:
                # If not then make an arbitrary choice from those available and increment the parsimony score
                node.final_tissue = random.choice(list(node.tissue_set))
                node.parsimony_score += 1
            node.name = f"{node.name}_{node.final_tissue}"
            # Recursively process children
            for child in node.children:
                preorder(child, node.final_tissue)
        else:
            # Check if leaf nodes are different tissues than their parents
            if parent_tissue != node.final_tissue:
                node.parsimony_score += 1

    # copy the input tree to avoid changing it in place
    copy_tree = tree.copy()

    # Run the postorder to get candidate tissues at each node
    postorder(copy_tree)

    # Initialize the parsimony scores
    for node in copy_tree.traverse():
        node.parsimony_score = 0

    # Assign the ancestral tissues for each node and update the parsimony score
    preorder(copy_tree)

    # Obtain the total parsimony score for the tree
    total_parsimony_score = sum(node.parsimony_score for node in copy_tree.traverse())

    return copy_tree

# User inputs
tree_file = sys.argv[1] # newick file
leaf_tissues_tsv = sys.argv[2]  # tsv file with tip cell names and tissue labels as columns, and no header
outdir = sys.argv[3] # where to write the output
primary_tissue = sys.argv[4] # the known tissue label of the root node

tree = ete3.Tree(tree_file, format=8)
tissue_map = pd.read_csv(leaf_tissues_tsv, sep=r'\s+', header=None, names=['cell', 'tissue'], dtype={'cell': str, 'tissue': str})

# get results
parsimony_tree = label_tissues_parsimony(tree, tissue_map)

parsimony_output = outdir + "/parsimony_tissues.nwk"
parsimony_tree.write(outfile = parsimony_output, format=8, format_root_node=True)
        