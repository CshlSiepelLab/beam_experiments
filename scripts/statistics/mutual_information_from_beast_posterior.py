#!/usr/bin/env python3

import dendropy
from collections import defaultdict
import numpy as np
import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor

def normalized_mutual_information(count_matrix):
    P = count_matrix / np.sum(count_matrix)
    p_x = np.sum(P, axis=1)
    p_y = np.sum(P, axis=0)
    
    # Compute mutual information
    MI = 0
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            if P[i,j] > 0:
                MI += P[i,j] * np.log2(P[i,j] / (p_x[i] * p_y[j]))

    # Compute entropies
    H_x = -np.sum(p_x[p_x > 0] * np.log2(p_x[p_x > 0]))
    H_y = -np.sum(p_y[p_y > 0] * np.log2(p_y[p_y > 0]))
    
    # Normalize by minimum entropy
    NMI = MI / min(H_x, H_y)
    
    return NMI

def process_tree(tree, origin_tissue):
    tissue_types = set()
    migration_counts = defaultdict(lambda: defaultdict(int))
    
    for node in tree.preorder_node_iter():
        if node.parent_node:
            source = node.parent_node.annotations.get_value('location')
            target = node.annotations.get_value('location')
            tissue_types.add(source)
            tissue_types.add(target)
            migration_counts[source][target] += 1
        else:
            target = node.annotations.get_value('location')
            tissue_types.add(origin_tissue)
            tissue_types.add(target)
            migration_counts[origin_tissue][target] += 1
    
    return tissue_types, migration_counts

def compute_migration_mutual_info(nexus_path, origin_tissue, threads):
    trees = dendropy.TreeList.get(path=nexus_path, schema='nexus')
    burnin_percent = 0.1
    num_beast_trees = len(trees)
    num_discard = round(num_beast_trees * burnin_percent)
    trees = trees[num_discard:]
    num_beast_trees = len(trees)
    
    tissue_types = set()
    migration_counts = defaultdict(lambda: defaultdict(int))
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(process_tree, trees, [origin_tissue] * len(trees))
    
    for t_types, m_counts in results:
        tissue_types.update(t_types)
        for source, targets in m_counts.items():
            for target, count in targets.items():
                migration_counts[source][target] += count
    
    tissue_list = sorted(tissue_types - {origin_tissue})
    tissue_list.insert(0, origin_tissue)
    n = len(tissue_list)
    count_matrix = np.zeros((n, n))
    for i, source in enumerate(tissue_list):
        for j, target in enumerate(tissue_list):
            count_matrix[i, j] = migration_counts[source][target]
    
    mutual_info = normalized_mutual_information(count_matrix)
    
    return mutual_info, count_matrix, tissue_list

# posterior_file = sys.argv[1]
# origin_tissue = sys.argv[2]
# outdir = sys.argv[3]
# threads = int(sys.argv[4])

# testing
posterior_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_1_20_25_uniform_50cells_50sites_data_7_24_24/beam_gtr/mS_854/combined.trees"
origin_tissue = "P"
outdir = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_1_20_25_uniform_50cells_50sites_data_7_24_24/beam_gtr/mS_854"
threads = 50

mi, counts, tissues = compute_migration_mutual_info(posterior_file, origin_tissue, threads)

with open(os.path.join(outdir, "posterior_trees_migration_count_matrix_for_mutual_information.csv"), mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow([""] + tissues)
    for i, row in enumerate(counts):
        writer.writerow([tissues[i]] + [int(value) for value in row])

with open(os.path.join(outdir, "posterior_trees_migration_mutual_information.txt"), mode='w') as file:
    file.write(str(mi))
