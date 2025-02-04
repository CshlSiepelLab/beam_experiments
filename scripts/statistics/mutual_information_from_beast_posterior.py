#!/usr/bin/env python3

import dendropy
from collections import defaultdict
import numpy as np
from sklearn.metrics import mutual_info_score

def compute_migration_mutual_info(nexus_path, origin_tissue):
    # Load trees from nexus file
    trees = dendropy.TreeList.get(path=nexus_path, schema='nexus')
    
    # Initialize count matrix 
    tissue_types = set()
    migration_counts = defaultdict(lambda: defaultdict(int))
    
    # Count migrations across all trees
    for tree in trees:
        for node in tree.preorder_node_iter():
            # for all nodes except the root, use the annotations to get the migration event
            if node.parent_node:
                source = node.parent_node.annotations.get_value('tissue')
                target = node.annotations.get_value('tissue') 
                tissue_types.add(source)
                tissue_types.add(target)
                migration_counts[source][target] += 1
            # for the root, use the input origin primary tissue the origin above the root which is not in the newick
            else:
                target = node.annotations.get_value('tissue')
                tissue_types.add(origin_tissue)
                tissue_types.add(target)
                migration_counts[origin_tissue][target] += 1

    
    # Convert to numpy array
    tissue_list = sorted(tissue_types)
    n = len(tissue_list)
    count_matrix = np.zeros((n,n))
    for i, source in enumerate(tissue_list):
        for j, target in enumerate(tissue_list):
            count_matrix[i,j] = migration_counts[source][target]
    
    # Calculate mutual information
    source_dist = count_matrix.sum(axis=1) / count_matrix.sum()
    target_dist = count_matrix.sum(axis=0) / count_matrix.sum()
    mutual_info = mutual_info_score(source_dist, target_dist)
    
    return mutual_info, count_matrix, tissue_list

posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_1_20_25_uniform_50cells_50sites_data_7_24_24/beam_gtr/mS_854/combined.trees"
origin_tissue = "P"

mi, counts, tissues = compute_migration_mutual_info(posterior_file, origin_tissue)

print(f"Mutual Information: {mi}")