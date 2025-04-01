#!/usr/bin/env python3

import os, sys
from ete3 import Tree

edge_list_file = sys.argv[1]

outfile = (
    os.path.dirname(edge_list_file)
    + "/"
    + os.path.basename(edge_list_file).split(".")[0]
    + ".nwk"
)

edges = []
with open(edge_list_file, "r") as file:
    for line in file:
        fields = line.strip().split()
        edge = (fields[0], fields[1])
        edges.append(edge)

tree = Tree.from_parent_child_table(edges)

print(outfile)
tree.write(outfile=outfile, format=8)
