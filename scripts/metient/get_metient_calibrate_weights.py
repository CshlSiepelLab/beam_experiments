
import os, sys
import gzip
import pickle
from metient import metient as met
import pandas as pd
import numpy as np

# Default weights that have been calibrated to real data
weights = met.Weights()
print("Default weights (not used, but printed to compare to here):")
print(weights.mig, weights.comig, weights.seed_site)

# Get new calibrated weights from simulated data
input_dir = "/mnt/stored_results/beam/latest_results/snakemake_performance_uniform_50cells_50sites_data_7_24_24/metient"
dataset_names = [os.path.basename(d) for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
print(len(dataset_names))

clone_tree_fns = [os.path.join(input_dir, name, "tree.txt") for name in dataset_names]
ref_var_fns = [os.path.join(input_dir, name, "metadata.tsv") for name in dataset_names]

print_config = met.PrintConfig(visualize=True, verbose=False, k_best_trees=5)

output_dir = "/mnt/stored_results/beam/latest_results/snakemake_performance_uniform_50cells_50sites_data_7_24_24/metient_calibrate_80_ideal_sims_8_12_25"
os.makedirs(output_dir, exist_ok=True)

met.calibrate_label_clone_tree(clone_tree_fns, ref_var_fns, print_config, output_dir, dataset_names, solve_polytomies=True)


