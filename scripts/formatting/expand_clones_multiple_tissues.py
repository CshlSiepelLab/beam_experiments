#!/usr/bin/env python3

import sys
import os
import pandas as pd

# inputs
matrix_file = sys.argv[1]
tissues_file = sys.argv[2]
output_dir = sys.argv[3]

# # testing
# matrix_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/billy_bladder_cancer_data_9_1_24/metastabayes/MMUS1782/CP00/temp_matrix.tsv"
# tissues_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/billy_bladder_cancer_data_9_1_24/metastabayes/MMUS1782/CP00/temp_tissues.tsv"
# output_dir ="./"

# read in matrix and tissues files to df where the files are tsv files
matrix_df = pd.read_csv(matrix_file, sep='\t', index_col=0)
tissues_df = pd.read_csv(tissues_file, sep='\t')

# only keep relevant columns group_name and tissues
tissues_df = tissues_df[['group_name', 'tissues']]

# find clones with more than one tissue
clones_with_multiple_tissues = tissues_df.loc[tissues_df['tissues'].str.contains(',', na=False), 'group_name'].values.tolist()

new_matrix_rows = []
new_tissues_rows = []
clones_to_drop = []

for clone in clones_with_multiple_tissues:
    tissues = tissues_df.loc[tissues_df['group_name'] == clone, 'tissues'].values[0].split(',')
    for i, tissue in enumerate(tissues):
        new_row = matrix_df.loc[matrix_df.index == clone].values[0]
        new_matrix_rows.append(pd.Series(new_row, index=matrix_df.columns, name=f'{clone}_{i}'))
        new_tissues_rows.append({'group_name': f'{clone}_{i}', 'tissues': tissue})
    clones_to_drop.append(clone)
    tissues_df = tissues_df.loc[tissues_df['group_name'] != clone]

# remove the original entries for the clone
for clone in clones_to_drop:
    matrix_df = matrix_df.drop(index=clone)

# concatenate new rows to the original DataFrames
matrix_df = pd.concat([matrix_df, pd.DataFrame(new_matrix_rows)])
tissues_df = pd.concat([tissues_df, pd.DataFrame(new_tissues_rows)], ignore_index=True)

# write outputs
matrix_df.to_csv(os.path.join(output_dir, "expanded_clones_matrix.tsv"), sep='\t')
tissues_df.to_csv(os.path.join(output_dir, "expanded_clones_tissues.tsv"), sep='\t', index=False)