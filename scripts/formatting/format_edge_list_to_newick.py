#!/usr/bin/env python3

import os, sys
from ete3 import Tree

def main():
    edge_list_file = sys.argv[1]

    outfile = edge_list_file.split(".")[0] + "edges_to_newick.nwk"

    edges=[]
    with open(edge_list_file, "r") as file:
        for line in file:
            fields = line.strip().split(" ")
            edge = set(fields[0], fields[1])
            edges.append(edge)
    
    tree = Tree.from_parent_child_table(edges)

    tree.write(outfile, format=8)

if __name__ == "__main__":
    main()