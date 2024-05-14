#!/usr/env/bin python3

import sys
# import math
import pandas as pd
import scipy.stats
import numpy as np

# # user input
# file_path = sys.argv[1]

# testing
file_path = 'results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/1983/1983_indel_character_matrix.tsv'

# read in mutation matrix
matrix = pd.read_csv(file_path, sep='\t', header=0, index_col=0)

# count the unique mutations and their frequencies per column
pis = 0
for col in matrix.columns:
    muts = matrix[col].values
    mut_counts = {}
    for mut in muts:
        if mut == 0:
            continue
        elif mut not in mut_counts:
            mut_counts[mut] = 1
        else:
            mut_counts[mut] += 1
    i = 0
    for mut, count in mut_counts.items():
        if count >= 2:
            i +=1
        if i >= 2:
            pis += 1
            break

print(f"Parsimony informative sites: {pis}")
