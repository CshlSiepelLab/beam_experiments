
import sys
import os

graph_file = sys.argv[1]    # Should be the G*parant_child.txt file
tissue_key_file = sys.argv[2]      # Should be the G*labels.txt file
primary_tissue = sys.argv[3]  # Should be the primary tissue name matching the label in the tissue_key_file, e.g. 'PRL' for Serio et al. data

# Get the mapping of tissue indices to tissue names since the machina output graph edges use indices
tissue_key = {}
with open(tissue_key_file, 'r') as f:
    for line in f.readlines():
        line = line.strip()
        if not line:
            continue
        index, tissue_name = line.split(' ')
        tissue_key[index] = tissue_name

# Read in edges from the graph file and check for PR (primary reseeding) or M2M (met to met) as a boolean for the full graph
m2m = False
pr = False
with open(graph_file, 'r') as f:
    for edge in f.readlines():
        edge = edge.strip()
        if not edge:
            continue
        source, target = edge.split(' ')
        # Convert source and target to tissue names using the tissue_key
        source_name = tissue_key[source]
        target_name = tissue_key[target]
        # Check for pr
        if source_name != primary_tissue and target_name == primary_tissue:
            pr = True
        # Check for m2m
        if source_name != primary_tissue and target_name != primary_tissue:
            m2m = True
        
        # If both pr and m2m are found, we can stop checking further
        if pr and m2m:
            break

# Print the results
print(f"pr: {pr}\nm2m: {m2m}")
