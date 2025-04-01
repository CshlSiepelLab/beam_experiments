#!/usr/bin/env python3

import sys
from ete3 import Tree
from dendropy import TreeList


def compute_rf(reference_tree, target_tree):
    # Helper function to compute RF distance between two trees
    def get_splits(tree):
        splits = set()
        for node in tree.traverse("postorder"):
            if not node.is_leaf():
                split = frozenset(leaf.name for leaf in node.get_leaves())
                splits.add(split)
        return splits

    ref_splits = get_splits(reference_tree)
    target_splits = get_splits(target_tree)

    # RF distance is the number of splits unique to each tree
    unique_to_ref = len(ref_splits - target_splits)
    unique_to_target = len(target_splits - ref_splits)
    total_splits = len(ref_splits) + len(target_splits)

    # Normalize RF distance by dividing by the total number of splits
    normalized_rf = (
        (unique_to_ref + unique_to_target) / total_splits if total_splits > 0 else 0
    )
    return normalized_rf


def compute_rf_distances(ref_tree, beam_trees):
    # Load the reference tree
    reference_tree = Tree(ref_tree, format=1)

    # Load the trees using DendroPy
    tree_list = TreeList.get(path=beam_trees, schema="nexus")

    # Convert the trees to Newick strings and load them into ete3 Tree objects
    posterior_trees = [
        Tree(tree.as_string(schema="newick"), format=1) for tree in tree_list
    ]

    rf_distances = []
    for tree in posterior_trees:
        rf = compute_rf(reference_tree, tree)
        rf_distances.append(rf)

    # Calculate the expectation of RF distances
    expectation_rf = sum(rf_distances) / len(rf_distances) if rf_distances else 0
    return expectation_rf, rf_distances


def main():
    ref_tree = sys.argv[1]
    beam_trees = sys.argv[2]
    outfile = sys.argv[3]

    expectation_rf, rf_distances = compute_rf_distances(ref_tree, beam_trees)
    with open(outfile, "w") as out:
        out.write(str(expectation_rf))


if __name__ == "__main__":
    main()
