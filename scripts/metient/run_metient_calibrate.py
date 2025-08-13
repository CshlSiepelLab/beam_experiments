
import os, sys
import gzip
import pickle
from metient import metient as met
import pandas as pd
import numpy as np
import random

# Default weights that have been calibrated to real data
weights = met.Weights()
print("Default weights (not used, but printed to compare to here):")
print(weights.mig, weights.comig, weights.seed_site)

# Get new calibrated weights from simulated data
input_dir = "/mnt/stored_results/beam/latest_results/snakemake_performance_uniform_50cells_50sites_data_7_24_24/metient"
dataset_names = [os.path.basename(d) for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
print(f"Found {len(dataset_names)} datasets:", dataset_names)

clone_tree_fns = [os.path.join(input_dir, name, "tree.txt") for name in dataset_names]
ref_var_fns = [os.path.join(input_dir, name, "metadata.tsv") for name in dataset_names]

print_config = met.PrintConfig(visualize=True, verbose=False, k_best_trees=5)

output_dir = "/mnt/stored_results/beam/latest_results/snakemake_performance_uniform_50cells_50sites_data_7_24_24/metient_calibrate_80_ideal_sims_8_12_25"
os.makedirs(output_dir, exist_ok=True)

# Run metient calibrate
met.calibrate_label_clone_tree(clone_tree_fns, ref_var_fns, print_config, output_dir, dataset_names, solve_polytomies=True)

# Get migration graphs for each dataset from metient results pkl file
for name in dataset_names:
    # Sort through results to obtain samples from the solution space
    with gzip.open(os.path.join(output_dir, "calibrate", f"{name}_P.pkl.gz"), "rb") as f:
        pckl = pickle.load(f)

    # Obtain samples from the solution space
    num_samples = 1024
    num_results = len(pckl["clone_tree_labelings"])
    if num_results < num_samples:
        num_samples = num_results

    tissues = pckl["ordered_anatomical_sites"]
    migration_graphs = np.empty(num_samples, dtype=dict)
    losses = np.empty(num_samples)

    for i in range(num_samples):
        V = pckl["clone_tree_labelings"][i]
        A = pckl["full_adjacency_matrices"][i]
        G = met.migration_graph(V, A)
        df_G = pd.DataFrame(G.numpy(), columns=tissues, index=tissues)
        df_dict = df_G.transpose().to_dict()
        migration_graphs[i] = df_dict
        losses[i] = pckl["losses"][i]

    with open(f"{output_dir}/calibrate/{name}_migration_graphs.txt", "w") as file:
        file.write(f"loss\tmigration_graph\n")
        for l, g in zip(losses, migration_graphs):
            file.write(f"{l}\t{g}\n")


