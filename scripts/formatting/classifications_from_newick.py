#!/usr/bin/env python 3

import sys
from ete3 import Tree

newick_file = sys.argv[1]
primary_tissue = sys.argv[2]
output_file = sys.argv[3]

tree = Tree(newick_file, format=8)

met_to_met = False
primary_reseeding = False

for node in tree.traverse():
    if node.is_root():
        continue
    else:
        node_tissue = node.name.split("_")[-1]
        parent_tissue = node.up.name.split("_")[-1]
        if node_tissue == primary_tissue and parent_tissue != primary_tissue:
            primary_reseeding = True
        if node_tissue != parent_tissue and node_tissue != primary_tissue and parent_tissue != primary_tissue:
            met_to_met = True

with open(output_file, "w") as f:
    f.write(f"met_to_met,primary_reseeding\n")
    f.write(f"{met_to_met},{primary_reseeding}\n")
