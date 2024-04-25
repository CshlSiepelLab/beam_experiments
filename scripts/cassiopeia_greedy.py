#!/usr/bin/env python3
import os, sys
import pandas as pd
import cassiopeia as cas
from ete3 import Tree


def main():
    character_matrix_tsv = sys.argv[1]

    # get dir path to input file for output to the same dir
    outprefix = os.path.dirname(character_matrix_tsv)

    # read in final matrix
    final_matrix = pd.read_csv(character_matrix_tsv, sep='\t', index_col=0)
    final_matrix.index = final_matrix.index.astype(str)

    # solve cassiopeia greedy
    reconstructed_tree = cas.data.CassiopeiaTree(character_matrix = final_matrix, missing_state_indicator = -1)
    greedy_solver = cas.solver.VanillaGreedySolver()
    greedy_solver.solve(reconstructed_tree)

    # make ete3 tree to write newick with internal node labels
    connections = reconstructed_tree.edges
    tree = Tree.from_parent_child_table(connections)

    # rename internal nodes from cassiopeia defaults
    i = 0
    for node in tree.traverse():
        if node.is_leaf() == True:
            continue
        else:
            node.name = f"node{i}"
            i+=1

    out_tree_infer = outprefix + "/cassiopeia_greedy_inferred.nwk"
    with open(out_tree_infer,'w') as it:
        it.write(tree.write(format=8))

if __name__ == "__main__":
    main()