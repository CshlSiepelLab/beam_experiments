#!/usr/bin/env python3

from ete3 import Tree
import sys


def classify_is_reseeding_from_newick(newick_file, primary_tissue, outfile):
    tree = Tree(newick_file, format=8)

    # make sure the root node is labeled with the primary tissue
    root = tree.get_tree_root()
    if not root.name.split("_")[-1] == primary_tissue:
        root.name = f"{tree.name}_{primary_tissue}"

    def traverse_and_check(node, primary_tissue):
        if node.is_root():
            return False
        parent_tissue = node.up.name.split("_")[-1]
        node_tissue = node.name.split("_")[-1]
        if node_tissue == primary_tissue and parent_tissue != primary_tissue:
            return True
        for child in node.children:
            if traverse_and_check(child, primary_tissue):
                return True
        return False

    result = "no"
    for node in tree.traverse():
        if traverse_and_check(node, primary_tissue):
            result = "yes"
            break

    with open(outfile, "w") as f:
        f.write(result)


newick_file = sys.argv[1]
primary_tissue = sys.argv[2]
outfile = sys.argv[3]

classify_is_reseeding_from_newick(newick_file, primary_tissue, outfile)
