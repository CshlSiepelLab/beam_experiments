#!/usr/bin/env python3

from ete3 import Tree
import sys
import ast

# inputs
metient = sys.argv[1]
primary_tissue = sys.argv[2]
outfile = sys.argv[3]

# get metient counts and boolean topology classifications
metient_loss_graphs = {}
with open(metient, "r") as f:
    next(f)  # Skip the first line header
    i=1
    for line in f.readlines():
        loss, graph = line.strip().split("\t")
        metient_loss_graphs[i] = graph
        i+=1

reseeding_counts = []

for loss, graph in metient_loss_graphs.items():
    graph_dict = ast.literal_eval(graph)

    reseeding = False

    for key, value in graph_dict.items():
        if value[primary_tissue] != 0:
            reseeding = True
            break

    reseeding_counts.append(reseeding)

metient_reseeding = max(set(reseeding_counts), key=reseeding_counts.count)

with open(outfile, "w") as f:
    if metient_reseeding:
        f.write("yes")
    else:
        f.write("no")