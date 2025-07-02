#!/usr/bin/env python3

import sys
import os
import random
import dendropy

# nexus_file = "/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25/beam_gtr/5k/26/combined_subsample.trees"

nexus_file = sys.argv[1]

outdir = os.path.dirname(nexus_file)
num_samples = 5

# Read trees from nexus file
tree_collection = dendropy.TreeList.get_from_path(nexus_file, schema="nexus")
trees = list(tree_collection)

# Discard first 50% of trees
half_point = len(trees) // 2
remaining_trees = trees[half_point:]

# Sample trees from the remaining 50%
sampled_trees = random.sample(remaining_trees, num_samples)

# Process each sampled tree to get a newick with tissue annotations appended to node names
for i, tree in enumerate(sampled_trees, 1):
    j = 0
    for node in tree.nodes():
        if node.is_leaf():
            continue
        else:
            tissue_label = str(node.annotations['location'].value).split("location=")[-1].replace("'", "")
            # Set node name to tissue.currentname format
            current_name = node.label
            if current_name != None:
                new_name = f"{tissue_label}.{current_name}"
                node.label = new_name
            else:
                new_name = f"{tissue_label}.node_{j}"
                node.label = new_name
                j += 1
            
    # Write to newick file
    outfile = os.path.join(outdir, f"sampled_tree_{i}_with_tissue_appended_to_name.nwk")
    tree.write(
        path=outfile,
        schema="newick",
        suppress_internal_node_labels=False
    )
