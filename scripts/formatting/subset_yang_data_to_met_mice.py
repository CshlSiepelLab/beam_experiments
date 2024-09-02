#!/usr/bin/env python3

import pandas as pd
import re
import cassiopeia as cas

###################
### PREPROCESSING
###################
allele_filepath = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/metadata/KPTracer.alleleTable.FINAL.txt"
outdir = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata"

allele_df = pd.read_csv(allele_filepath, sep='\t', index_col = 0)


###################
### FIND MET MICE
###################

# make new column with tissue labels only
allele_df['tissue'] = [''.join(filter(str.isalpha, name[2])) for name in allele_df['Tumor'].str.split('_')]

# get met lineage names and tissues
tumor_names = allele_df.groupby('MetFamily')['tissue'].unique()

# keep only lineages with both primary and met tissues
tumor_names = tumor_names[tumor_names.apply(lambda x: any(name.startswith("T") for name in x) and any(not name.startswith("T") for name in x))]

# keep only subset MetFamily mice
allele_df = allele_df[allele_df['MetFamily'].isin(tumor_names.index)]

###################
### MAKE CHAR MATRIX
###################

# write mouse specific character matrix to file
for MetFamily, df in allele_df.groupby('MetFamily'):

    # calculate the threshold to make sure sites are only included if at least one cell has a different mutation than the rest, i.e. remove uninformative sites with 100% of the site as the same mutation
    cell_threshold = 1 - (1 / df['cellBC'].nunique())

    # built in cassiopeia function to convert allele table to character matrix
    char_matrix = cas.pp.convert_alleletable_to_character_matrix(df, allele_rep_thresh = cell_threshold)

    char_matrix_df = char_matrix[0]
    mut_dict = char_matrix[2]

    # convert char matrix to successive char matrix
    successive_char_matrix = char_matrix_df.copy()
    successive_mut_dict = {}
    i = 1
    for clone, row in char_matrix_df.iterrows():
        for site, mut in row.items():
            # skip unedited or missing sites
            mut = int(mut)
            if mut == 0 or mut == -1:
                continue
            mut_str = mut_dict[int(site[1:])-1][mut]
            # replace the mutation with the successive mutation
            if mut_str not in successive_mut_dict:
                successive_mut_dict[mut_str] = i
                new_mut_value = i
                i += 1
            else:
                new_mut_value = successive_mut_dict[mut_str]
            successive_char_matrix.loc[clone, site] = new_mut_value

    # rename columns of successive char matrix to be successive themselves
    successive_char_matrix.columns = [f"r{i}" for i in range(1, len(successive_char_matrix.columns)+1)]

    # output successive char matrix for all cells
    mouse_outfile = f"{outdir}/{MetFamily}_successive_char_matrix.txt"
    successive_char_matrix.index.name = 'cellBC'
    # successive_char_matrix.to_csv(mouse_outfile, sep="\t", index=True)

    # write mutation dictionary to file
    mut_dict_outfile = mouse_outfile.replace(".txt", f"_mut_dict.txt")
    # with open(mut_dict_outfile, "w") as f:
    #     f.write(f"mut_id\tmut_str\n")
    #     for str, id in successive_mut_dict.items():
    #         f.write(f"{id}\t{str}\n")

    # collapse the cells to only unique rows and output collapsing dict of cellBCs and tissue labels
    all_columns = successive_char_matrix.columns.tolist()
    sorted_char_matrix = successive_char_matrix.sort_values(by=all_columns)
    unique_rows = sorted_char_matrix.drop_duplicates(keep='first')
    group_names = [f"clone{i+1}" for i in range(len(unique_rows))]
    group_to_originals = {}
    group_to_tissues = {}
    for group_name, (_, unique_row) in zip(group_names, unique_rows.iterrows()):
        # Find all rows in sorted_char_matrix that match the unique_row
        original_row_names = sorted_char_matrix[sorted_char_matrix.eq(unique_row).all(axis=1)].index.tolist()
        group_to_originals[group_name] = original_row_names
        original_tissues = set(df[df['cellBC'].isin(original_row_names)]['tissue'].values.tolist())
        group_to_tissues[group_name] = original_tissues
    
    # replace index names in unique_rows with the appropriate group name
    unique_rows.index = group_names

    print(MetFamily, f"cells: {len(successive_char_matrix)}", f"clones: {len(unique_rows)}", f"sites: {len(successive_char_matrix.columns)}")

    # write unique rows to file
    unique_rows_outfile = mouse_outfile.replace(".txt", f"_collapsed.txt")
    # unique_rows.to_csv(unique_rows_outfile, sep="\t", index=True)

    # write collapsing dict of cellBCs and tissue labels to file
    collapsing_dict_outfile = mouse_outfile.replace(".txt", f"_collapsing_dict.txt")
    # with open(collapsing_dict_outfile, "w") as f:
    #     f.write(f"group_name\tcellBCs\ttissues\n")
    #     for group_name in group_names:
    #         cellBCs = ','.join(list(group_to_originals[group_name]))
    #         tissues = ','.join(list(group_to_tissues[group_name]))
    #         f.write(f"{group_name}\t{cellBCs}\t{tissues}\n")
