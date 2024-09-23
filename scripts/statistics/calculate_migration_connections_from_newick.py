#!/usr/bin/env python3

import sys
from ete3 import Tree


newick_file = sys.argv[1]
# newick_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_9_12_24/raw_data/817/tissue_labeled_tree.nwk"
primary_tissue = "P"

tree = Tree(newick_file, format=1)

# set origin to known primary
tree.name = f"root_{primary_tissue}"

connections = 0

connections_done = []

for node in tree.traverse():
    if node.is_leaf():
        continue
    else:
        name, tissue = node.name.split("_")
        for child in node.get_children():
            child_name, child_tissue = child.name.split("_")
            if tissue != child_tissue and f"{name}_{child_name}" not in connections_done:
                connections_done.append(f"{name}_{child_name}")
                connections += 1

print(f"{connections}")

