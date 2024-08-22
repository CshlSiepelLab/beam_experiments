#!/usr/bin/env python3

import sys
import os
import pandas as pd


# # user inputs
# cas_char_matrix_files = sys.argv[1]
# cas_mut_dict_files = sys.argv[2]

# make sure the ordering is matched between char matrix and mut dict input files
cas_char_matrix_files = ["/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3513_NT_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3466_Lkb1_T2_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3515_Lkb1_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3508_Apc_T2_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3451_Lkb1_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3457_Apc_T4_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3724_NT_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3454_Lkb1_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3465_Lkb1_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3457_Apc_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3460_Lkb1_T1_char_matrix_collapsed.txt",
                        "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3519_Lkb1_T1_char_matrix_collapsed.txt"]

cas_mut_dict_files = ["/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3513_NT_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3466_Lkb1_T2_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3515_Lkb1_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3508_Apc_T2_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3451_Lkb1_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3457_Apc_T4_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3724_NT_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3454_Lkb1_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3465_Lkb1_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3457_Apc_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3460_Lkb1_T1_char_matrix_mut_dict.txt",
                    "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered_3519_Lkb1_T1_char_matrix_mut_dict.txt"]

for cas_char_matrix, cas_mut_dict in zip(cas_char_matrix_files, cas_mut_dict_files):

    cas_char_matrix_df = pd.read_csv(cas_char_matrix, sep='\t', index_col=0)
    cas_mut_dict_df = pd.read_csv(cas_mut_dict, sep='\t', index_col=0)

    successive_char_matrix = cas_char_matrix_df.copy()
    successive_mut_dict = {}
    i = 1
    for clone, row in cas_char_matrix_df.iterrows():
        for site, mut in row.items():
            # skip unedited or missing sites
            if mut == 0 or mut == -1:
                continue
            subset_mut_dict_df = cas_mut_dict_df.loc[[site], :]
            mut_str = subset_mut_dict_df.loc[subset_mut_dict_df['mut_id'] == mut, 'mut_str'].values[0]
            # replace the mutation with the successive mutation
            if mut_str not in successive_mut_dict:
                successive_mut_dict[mut_str] = i
                new_mut_value = i
                i += 1
            else:
                new_mut_value = successive_mut_dict[mut_str]
            successive_char_matrix.loc[clone, site] = new_mut_value

    outfile_char_matrix = cas_char_matrix.replace(".txt", "_successive_across_sites.txt")
    outfile_mut_dict = cas_mut_dict.replace(".txt", "_successive_across_sites.txt")

    successive_char_matrix.to_csv(outfile_char_matrix, sep='\t')
    with open(outfile_mut_dict, 'w') as f:
        f.write("mut_id\tmut_str\n")
        for mut_str, mut_id in successive_mut_dict.items():
            f.write(f"{mut_id}\t{mut_str}\n")



