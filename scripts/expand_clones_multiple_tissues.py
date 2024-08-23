#!/usr/bin/env python3

import sys
import os
import pandas as pd

# inputs
matrix_file = sys.argv[1]
tissues_file = sys.argv[2]
output_dir = sys.argv[3]

# # testing
# matrix_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3451_Lkb1_T1_char_matrix_collapsed_successive_across_sites.txt"
# tissues_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3451_Lkb1_T1_char_matrix_collapsing_dict.txt"
# output_dir ="./"

# read in matrix and tissues files to df where the files are tsv files
matrix_df = pd.read_csv(matrix_file, sep='\t', index_col=0)
tissues_df = pd.read_csv(tissues_file, sep='\t')

# drop second column from tissues_df which contains the cell barcodes for each clone that were collapses into the clone
tissues_df = tissues_df.drop(columns=tissues_df.columns[1])

# find clones with more than one tissue
clones_with_multiple_tissues = tissues_df.loc[tissues_df['tissues'].str.contains(',', na=False), 'group_name'].values.tolist()

new_matrix_rows = []
new_tissues_rows = []

for clone in clones_with_multiple_tissues:
    tissues = tissues_df.loc[tissues_df['group_name'] == clone, 'tissues'].values[0].split(',')
    for i, tissue in enumerate(tissues):
        new_row = matrix_df.loc[matrix_df.index == clone].values[0]
        new_matrix_rows.append(pd.Series(new_row, index=matrix_df.columns, name=f'{clone}_{i}'))
        new_tissues_rows.append({'group_name': f'{clone}_{i}', 'tissues': tissue})
    # remove the original entries for the clone
    matrix_df = matrix_df.drop(index=clone)
    tissues_df = tissues_df.loc[tissues_df['group_name'] != clone]

# concatenate new rows to the original DataFrames
matrix_df = pd.concat([matrix_df, pd.DataFrame(new_matrix_rows)])
tissues_df = pd.concat([tissues_df, pd.DataFrame(new_tissues_rows)], ignore_index=True)

# write outputs
matrix_df.to_csv(os.path.join(output_dir, "expanded_clones_matrix.tsv"), sep='\t')
tissues_df.to_csv(os.path.join(output_dir, "expanded_clones_tissues.tsv"), sep='\t', index=False)