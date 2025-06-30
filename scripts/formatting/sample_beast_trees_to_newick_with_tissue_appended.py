#!/usr/bin/env python3

import sys
import os
import random
import dendropy

# Set parameters
nexus_file = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/26/combined.trees"
outdir = os.path.dirname(nexus_file)
num_samples = 10

# Read trees from nexus file
tree_collection = dendropy.TreeList.get_from_path(nexus_file, schema="nexus")
trees = list(tree_collection)

# Discard first 50% of trees
half_point = len(trees) // 2
remaining_trees = trees[half_point:]

# Sample trees from the remaining 50%
sampled_trees = random.sample(remaining_trees, num_samples)

# Process each sampled tree
for i, tree in enumerate(sampled_trees, 1):
    
    # Extract tissue label from all nodes in the tree
    j = 0
    for node in tree.nodes():
        # Check for tissue annotation in node annotations
        if 'location' in node.annotations:
            tissue_label = node.annotations['location']
            # Set node name to tissue.currentname format
            current_name = node.taxon.label if node.taxon else node.label
            if current_name:
                new_name = f"{tissue_label}.{current_name}"
                node.taxon.label = new_name
                node.label = new_name
            else:
                new_name = f"{tissue_label}.node_{j}"
                node.taxon.label = new_name
                node.label = new_name
                j += 1
    
    # Write to file
    outfile = os.path.join(outdir, f"sampled_tree_{i:03d}_with_tissue_appended_to_name.nwk")
    tree.write(
        path=outfile,
        schema="newick",
        suppress_internal_node_labels=False
    )
