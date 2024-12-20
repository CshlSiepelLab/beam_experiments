#!/usr/bin/env python3

import sys
import ast
import os
import ete3
from Bio import Phylo
from io import StringIO


def getMigrationComigrationCountsFromNewick(newick):
    print(newick)
    tree = ete3.Tree(newick, format=8)

    root = tree.get_tree_root()
    if not root.name:
        tree.get_tree_root().name = f"root_{primaryTissue}"

    migration_count = 0
    comigration_count = 0
    num_multiedges = 0
    met_to_met = False
    reseeding = False
    edges = []
    multiedges_already_checked = []

    for node in tree.traverse():
        if node.is_root():
            continue
        if node.up.is_root():
            parent_tissue = primaryTissue
        else:
            parent_tissue = node.up.name.split("_")[-1]
        child_tissue = node.name.split("_")[-1]
        if parent_tissue != child_tissue:
            migration_count += 1
            edge = f"{parent_tissue}_{child_tissue}"
            if edge not in edges: # All unique edges are considered a co-migration
                comigration_count += 1
            if edge in edges and edge not in multiedges_already_checked:
                num_multiedges += 1
                multiedges_already_checked.append(edge)
            if not met_to_met and parent_tissue != primaryTissue and child_tissue != primaryTissue:
                met_to_met = True
            if not reseeding and child_tissue == primaryTissue:
                reseeding = True
            edges.append(edge)
    
    if num_multiedges != 0:
        clonality = "Polyclonal"
    else:
        clonality = "Monoclonal"

    return migration_count, comigration_count, num_multiedges, met_to_met, reseeding, clonality


sim_name = sys.argv[1]
newick = sys.argv[2]
primaryTissue = sys.argv[3]
outfile = sys.argv[4]

migration_count, comigration_count, num_multiedges, met_to_met, reseeding, clonality = getMigrationComigrationCountsFromNewick(newick)

with open(outfile, "a") as f:
    f.write(f"{sim_name},{migration_count},{comigration_count},{num_multiedges},{met_to_met},{reseeding},{clonality}\n")

