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
            node.taxon.label = prediction
        except Exception as e:
            prediction = f"node{i}" + "_" + node.annotations.get_value('location')
            i += 1
        node.label = prediction
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("\'", ""), format=3)
    return ete_tree

def main():
    beast_file = sys.argv[1]

    #beast_file = "sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/joint_inference_beast_tissues.tree"

    out_file = beast_file + ".nwk"

    beast_tree = dendropy.Tree.get(path=beast_file, schema='nexus')
    remove_zero_length_nodes(beast_tree)
    beast_tree_ete = dendropy_beast_to_ete_newick_with_strict_locations(beast_tree)

    beast_tree_ete.write(outfile=out_file, format=8)

if __name__ == "__main__":
    main()
