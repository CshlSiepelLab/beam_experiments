#!/usr/bin/env python3

import sys
import os
import pandas as pd

char_matrix_file = sys.argv[1]
mut_dict_file = sys.argv[2]
outdir = sys.argv[3]

# # testing
# char_matrix_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/billy_bladder_cancer_data_9_1_24/metastabayes/MMUS1834/CP00/temp_matrix.tsv"
# mut_dict_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/billy_bladder_cancer_data_9_1_24/metastabayes/MMUS1834/mutation_dict.tsv"
# outdir = "./"

char_matrix_df = pd.read_csv(char_matrix_file, sep="\t", index_col=0)
mut_dict_df = pd.read_csv(mut_dict_file, sep="\t", index_col=0, header=None)

successive_char_matrix = char_matrix_df.copy()
successive_mut_dict = {}
i = 1
for clone, row in char_matrix_df.iterrows():
    for site, mut in row.items():
        mut = int(mut)
        # skip unedited or missing sites
        if mut == 0 or mut == -1:
            continue
        mut_str = mut_dict_df.loc[mut, 1]
        # replace the mutation with the successive mutation
        if mut_str not in successive_mut_dict:
            successive_mut_dict[mut_str] = i
            new_mut_value = i
            i += 1
        else:
            new_mut_value = successive_mut_dict[mut_str]
        successive_char_matrix.loc[clone, site] = new_mut_value

outfile_char_matrix = f"{outdir}/successive_char_matrix.csv"
outfile_mut_dict = f"{outdir}/successive_mut_dict.csv"

successive_char_matrix.to_csv(outfile_char_matrix, sep="\t")
with open(outfile_mut_dict, "w") as f:
    f.write("mut_id\tmut_str\n")
    for mut_str, mut_id in successive_mut_dict.items():
        f.write(f"{mut_id}\t{mut_str}\n")
