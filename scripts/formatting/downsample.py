#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
import random

def jaccard_distance(set1, set2):
    # Symmetric difference of two sets
    Symmetric_difference = set1.symmetric_difference(set2)

    # Unions of two sets
    union = set1.union(set2)
     
    return len(Symmetric_difference)/len(union)

def compute_distance(args):
    i, j = args
    return i, j, jaccard_distance(set(char_matrix.loc[i]), set(char_matrix.loc[j]))


# inputs
char_matrix_file = sys.argv[1]
tissue_labels_file = sys.argv[2]
threshold = float(sys.argv[3])
cores = int(sys.argv[4])
outprefix = sys.argv[5]

# # testing
# char_matrix_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/3451_Lkb1_T1_successive_char_matrix_collapsed.txt"
# tissue_labels_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/3451_Lkb1_T1_successive_char_matrix_collapsing_dict.txt"
# threshold = 0.1
# cores = 10
# outprefix = "test"

# read in data
char_matrix_original = pd.read_csv(char_matrix_file, sep='\t', index_col=0)
tissue_labels_original = pd.read_csv(tissue_labels_file, sep='\s+', index_col=0)

# remove any clones with more than one tissue label, which by default are informative and should be kept
tissue_labels = tissue_labels_original[tissue_labels_original['tissues'].str.contains(',', na=False) == False]
char_matrix = char_matrix_original.loc[tissue_labels.index]

# Do not include missing data
char_matrix = char_matrix.replace(-1, 0)

# Compute pairwise distances between rows
distance_matrix = pd.DataFrame(index=char_matrix.index, columns=char_matrix.index)

pool = Pool(processes=cores)
output = pool.map(compute_distance, [(i, j) for i in char_matrix.index for j in char_matrix.index])
pool.close()
pool.join()

for i, j, distance in output:
    distance_matrix.loc[i, j] = distance


# Find cells within threshold distance with the same tissue
tissues = list(set([t for tis in tissue_labels['tissues'] for t in tis.split(",")]))

clone_tissues = {}
for clone in tissue_labels.index:
    clone_tissues[clone] = tissue_labels.loc[clone, 'tissues'].split(",")

print("Initial number of clones: ", len(char_matrix.index))

filtered_clones = []
for tissue in tissues:
    tissue_clones = [clone for clone in clone_tissues if tissue in clone_tissues[clone]]
    for i in tissue_clones:
        if i in filtered_clones:
            continue
        for j in tissue_clones:
            if i == j:
                continue
            if i in filtered_clones or j in filtered_clones:
                continue

            if distance_matrix.loc[i, j] < threshold:
                remove = random.choice([i, j])
                filtered_clones.append(remove)
                print(f"Filtering out {remove} from {tissue} because of distance {distance_matrix.loc[i, j]} between {i} and {j}")

# Keep only clones not in filtered clones for both the char-matrix and tissue labels
char_matrix_filtered = char_matrix.loc[~char_matrix.index.isin(filtered_clones)]
tissue_labels_filtered = tissue_labels.loc[~tissue_labels.index.isin(filtered_clones)]
print("Filtered number of clones: ", len(char_matrix_filtered.index))

# write to files
char_matrix_filtered.to_csv(outprefix + "_char_matrix.tsv", sep='\t')
tissue_labels_filtered.to_csv(outprefix + "_tissue_labels.tsv", sep=',', header=False)

