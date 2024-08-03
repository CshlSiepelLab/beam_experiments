#!/usr/bin/env

import os, sys
import gzip
import pickle
from metient import metient as met
import pandas as pd
import numpy as np

tsv = sys.argv[1]
tree = sys.argv[2]
patient = sys.argv[3]
primary = sys.argv[4]
output_dir = sys.argv[5]

# tsv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/metient/2945/metadata.tsv"
# tree = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/metient/2945/tree.txt"
# patient = "2945"
# primary = "P"
# output_dir = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/metient/2945"

df = pd.read_csv(tsv, sep="\t")

print_config = met.PrintConfig(visualize=True, verbose=True, k_best_trees=5)
weights = met.Weights() # Use default weights which have been calibrated to real data

met.evaluate_label_clone_tree(tree, tsv, weights, print_config, output_dir, patient, solve_polytomies=True)


### Sort through results to obtain samples from the solution space
with gzip.open(os.path.join(output_dir, f"{patient}_{primary}.pkl.gz") ,"rb") as f:
    pckl = pickle.load(f)
# print(pckl.keys())

# obtain samples from the solution space up to 1024 samples if available
num_samples = 1024
num_results = len(pckl['clone_tree_labelings'])
if num_results < num_samples:
    num_samples = num_results

tissues = pckl['ordered_anatomical_sites']
migration_graphs = np.empty(num_samples, dtype=dict)
losses = np.empty(num_samples)
# trees = np.empty(num_samples,dtype=dict)


for i in range(num_samples):
    # obtain the migration graph
    V = pckl['clone_tree_labelings'][i]
    A = pckl['full_adjacency_matrices'][i]
    G = met.migration_graph(V, A)
    df_G = pd.DataFrame(G.numpy(), columns=tissues, index=tissues)
    df_dict = df_G.transpose().to_dict()
    migration_graphs[i] = df_dict
    losses[i] = pckl['losses'][i]

    # would be nice to record the full tree node labels here, but unsure how to do it given that some scenarios may have differing trees due to polytomy resolution, so will leave it blank to prevent pipeline bugs

# output samples to files
with open(f"{output_dir}/{patient}_{primary}_migration_graphs.txt", "w") as file:
    file.write(f"loss\tmigration_graph\n")
    for l, g in zip(losses, migration_graphs):
        file.write(f"{l}\t{g}\n")

# with open(f"{output_dir}/{patient}_{primary}_labeled_trees.txt", "w") as file:
#     for t in trees:
#         file.write(f"{t}\n")

    
    
