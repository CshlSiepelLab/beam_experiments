#!/usr/bin/env

import os, sys
from metient import metient as met
import pandas as pd

tsv = sys.argv[1]
tree = sys.argv[2]
patient = sys.argv[3]
output_dir = sys.argv[4]

df = pd.read_csv(tsv, sep="\t")

print_config = met.PrintConfig(visualize=True, verbose=False, k_best_trees=5)
weights = met.Weights(mig=1, comig=50, seed_site=100)

clone_tree_fn = os.path.join(tree)
ref_var_fn = os.path.join(tsv)
met.evaluate_label_clone_tree(clone_tree_fn, ref_var_fn, weights, print_config, 
                                output_dir, patient, solve_polytomies=True)

