#!/usr/bin/env python3

import sys
import pandas as pd
import cassiopeia as cas

def convert_matrix_to_successive(character_matrix, mut_dict, indel_priors):
    # Convert character matrix to successive character matrix
    successive_char_matrix = character_matrix.copy()
    successive_mut_dict = {}
    i = 1
    for clone, row in character_matrix.iterrows():
        for site, mut in row.items():
            mut = int(mut)

            # Skip undedited and missing sites
            if mut == 0 or mut == -1:
                continue

            mut_str = mut_dict[int(site[1:])-1][mut]
            # Replace the mutation with the successive mutation
            if mut_str not in successive_mut_dict:
                successive_mut_dict[mut_str] = i
                new_mut_value = i
                i += 1
            else:
                new_mut_value = successive_mut_dict[mut_str]
            successive_char_matrix.loc[clone, site] = new_mut_value
    
    # Compute edit rates for the successive matrix values
    successive_edit_rates = {}
    for mut_str, mut_int in successive_mut_dict.items():
        successive_edit_rates[mut_int] = indel_priors.loc[mut_str]['freq']

    return successive_char_matrix, successive_mut_dict, successive_edit_rates

infile = sys.argv[1]
lineage = int(sys.argv[2])
outdir = sys.argv[3]

# # testing
# infile = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/quinn_2021_real_data/GSE161363/GSM4905334_alleleTable.5k.txt'
# lineage=21
# outdir = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_1_22_25/successive_raw_data/5k'

# read in the proivded allele table
allele_table = pd.read_csv(infile, sep='\t', usecols = ['cellBC', 'intBC', 'r1', 'r2', 'r3', 'allele', 'LineageGroup', 'sampleID', 'readCount', 'UMI'])

# get indel priors as per Cassiopeis docs
indel_priors = cas.pp.compute_empirical_indel_priors(allele_table, grouping_variables=['intBC', 'LineageGroup'])

group = allele_table[allele_table['LineageGroup'] == lineage]

char_matrix_df, priors, mut_dict = cas.pp.convert_alleletable_to_character_matrix(group, missing_data_state='-1', allele_rep_thresh=0.95, mutation_priors = indel_priors)
successive_matrix, new_mut_dict, successive_edit_rates = convert_matrix_to_successive(char_matrix_df, mut_dict, indel_priors)

# write successive edit rates to a file
successive_edit_rates = dict(sorted(successive_edit_rates.items()))
with open(f"{outdir}/{lineage}_mutation_priors.txt", 'w') as f:
    f.write(f"mutation_code,rate\n")
    for key, value in successive_edit_rates.items():
        f.write(f"{key},{value}\n")

# write each lineage group's successive matrix to its own file
char_matrix_df.to_csv(f"{outdir}/{lineage}_original_character_matrix.tsv", sep='\t', index=True, header=True)

with open(f"{outdir}/{lineage}_original_chracter_int_to_mutation_string_dict.txt", 'w') as f:
    f.write(f"site_num,char_int,mut_str\n")
    for key, value in mut_dict.items():
        for k, v in value.items():
            f.write(f"{key},{k},{v}\n")

successive_matrix.to_csv(f"{outdir}/{lineage}_successive_character_matrix.tsv", sep='\t', index=True, header=True)

with open(f"{outdir}/{lineage}_successive_int_to_mutation_string_dict.txt", 'w') as f:
    f.write(f"successive_char_int,mut_str\n")
    for key, value in new_mut_dict.items():
        f.write(f"{value},{key}\n")

