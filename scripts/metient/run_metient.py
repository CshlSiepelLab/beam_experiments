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
weights = met.Weights() # Use default weights which have been calibrated to real data

met.evaluate_label_clone_tree(tree, tsv, weights, print_config, output_dir, patient, solve_polytomies=True)

