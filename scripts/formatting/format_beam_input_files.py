#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.Newick import Tree


indel_matrix_file=sys.argv[1]
outdir=sys.argv[2]

# assumes the output dir is named with the sim number
outname = os.path.basename(outdir)

# read in the indel matrix
indel_matrix = pd.read_csv(indel_matrix_file, sep="\t", index_col=0)

# replace mutation values with sequential values required for the editing model
done = []
mut_dict = {-1: -1} # keep as a dropout site to replace later
i = 1
for vals in indel_matrix.values.tolist():
    for v in vals:
        if v == -1 or v == 0 or v in mut_dict.keys():
            continue
        else:
            mut_dict[v] = i
            i = i + 1

# replace all entries in the indel_matrix with the mut_dict value
indel_matrix = indel_matrix.replace(mut_dict)

# get all indel proportions
muts = np.array([v for vals in indel_matrix.values.tolist() for v in vals if v != 0 and v != -1])
ordered_value_counts = np.unique(muts, return_counts=True)[1]
sum_proportions = sum(ordered_value_counts)
proportions = [str(count / sum_proportions) for count in ordered_value_counts]

# replace the -1 with the largest value + 1 for dropout as the last column in tidetree
max_val = max(mut_dict.values())
mut_dict[-1] = max_val + 1

# write mutation proportions for initial states in the edit model
outfile_proportions = f"{outdir}/{outname}_edit_rate_proportions.txt"
with open(outfile_proportions, "w") as file:
    file.write(" ".join(proportions))

# write mut dict to file
outfile_mut_dict = f"{outdir}/{outname}_original_mut_to_edit_model_mut.csv"
with open(outfile_mut_dict, "w") as file:
    file.write(f"original_mut,new_mut\n")
    for key,value in mut_dict.items():
        file.write(f"{key},{value}\n")

# write fasta file
outfile_fasta = f"{outdir}/{outname}.fasta"
with open(outfile_fasta, "w") as file:
    for index, row in indel_matrix.iterrows():
        sequence = ",".join(str(x) for x in row)
        file.write(f">{index}\n{sequence}\n")


