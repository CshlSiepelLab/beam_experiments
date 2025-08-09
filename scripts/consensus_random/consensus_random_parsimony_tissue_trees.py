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
            tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
            node.name = f"{node.name}_{tissue}"
        # internal nodes are chosen by consensus of leaf tissues
        else:
            node_name = node.name
            children = [
                tissues_df.loc[tissues_df["cell"] == str(leaf.name), "tissue"].values[0]
                for leaf in node.get_leaves()
            ]
            counts = Counter(children)
            most_common_elements = counts.most_common(2)
            # if there is a clear winner, choose that tissue
            if (
                len(most_common_elements) == 1
                or most_common_elements[0][1] != most_common_elements[1][1]
            ):
                consensus_tissue = most_common_elements[0][0]
            else:
                # tie breaker goes to primary
                if any(primary_tissue in elem for elem in most_common_elements):
                    consensus_tissue = primary_tissue
                # if tie doesn't involve primary, choose randomly
                else:
                    tied_elements = [elem[0] for elem in most_common_elements]
                    consensus_tissue = random.choice(tied_elements)
            # rename with consensus tissue label
            node.name = f"{node_name}_{consensus_tissue}"
    return copy_tree


def label_tissues_random(tree, tissues_df):
    copy_tree = tree.copy()
    tissues = tissues_df["tissue"].unique().tolist()
    for node in copy_tree.traverse():
        # leaves are assumed to be known labels
        if node.is_leaf():
            tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
            node.name = f"{node.name}_{tissue}"
        # internal nodes are chosen randomly
        else:
            random_tissue = random.choice(tissues)
            node.name = f"{node.name}_{random_tissue}"
    return copy_tree


def label_tissues_parsimony(tree, tissues_df):
    """
    Fitch parsimony algorithm to infer ancestral states of internal nodes for a tree
    """

    def postorder(node):
        if node.is_leaf():
            # Assign known tissue type from the tissues_df to the leaf node
            node.final_tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
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
                node.final_tissue = "P"
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
tree_file = sys.argv[1]
leaf_tissues_tsv = sys.argv[2]
outdir = sys.argv[3]

primary_tissue = "P"

tree = ete3.Tree(tree_file, format=8)

# Force the tree to be binary (ie. resovle polytomies)
tree.resolve_polytomy()

tissue_map = pd.read_csv(
    leaf_tissues_tsv,
    sep=r"\s+",
    header=None,
    names=["cell", "tissue"],
    dtype={"cell": str, "tissue": str},
)

# get results
random_tree = label_tissues_random(tree, tissue_map)
consensus_tree = label_tissues_consensus(tree, tissue_map)
parsimony_tree = label_tissues_parsimony(tree, tissue_map)

# output trees to newick files
random_output = outdir + "/random_tissues.nwk"
random_tree.write(outfile=random_output, format=8, format_root_node=True)

consensus_output = outdir + "/consensus_tissues.nwk"
consensus_tree.write(outfile=consensus_output, format=8, format_root_node=True)

parsimony_output = outdir + "/parsimony_tissues.nwk"
parsimony_tree.write(outfile=parsimony_output, format=8, format_root_node=True)
