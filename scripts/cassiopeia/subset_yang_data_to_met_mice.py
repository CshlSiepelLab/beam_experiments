#!/usr/bin/env python3

import pandas as pd
import re
import cassiopeia as cas

all_mice_metadata = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/metadata/KPTracer_meta.csv"
outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer_meta_FILTERED.csv"

metadata = pd.read_csv(all_mice_metadata, sep=",")

# drop rows with no tumor label read in as NaN
metadata.dropna(subset=['Tumor'], inplace=True)

# total number of mice in the raw data
original_total_mice = len(metadata['Mouse'].unique())
print(f"Total number of mice in the raw data: {original_total_mice}")

# total mice by genotype
original_total_mice_genotype = metadata.groupby('genotype')['Mouse'].nunique()
print(f"Total mice by genotype original: {original_total_mice_genotype}")

# make new column with only the anatomical site from the tumor label
metadata['tissue'] = [''.join(filter(str.isalpha, name)) for name in metadata['Tumor'].str.split('_').str[2]]

# get tumor names for each mouse
tumor_names = metadata.groupby('Mouse')['tissue'].unique()

# remove "Normal" entries in each group of tumor names
tumor_names = tumor_names.apply(lambda x: [name for name in x if name != "Normal"])

# subset tumor names to keep only groups with "T" and non-"T" entries
subset_tumor_names = tumor_names[tumor_names.apply(lambda x: "T" in x and sum(name != "T" for name in x) >= 2)]

# filter metadata based on subset tumor names
subset_metadata = metadata[metadata['Mouse'].isin(subset_tumor_names.index)]

# total number of mice in the filtered data
filtered_total_mice = len(subset_metadata['Mouse'].unique())
print(f"Total number of mice in the filtered data: {filtered_total_mice}")

# total mice by genotype
filtered_total_mice_genotype = subset_metadata.groupby('genotype')['Mouse'].nunique()
print(f"Total mice by genotype filtered: {filtered_total_mice_genotype}")

# get unique tissues
unique_tissues = set([item for sublist in subset_tumor_names.values.tolist() for item in sublist])

# remove any rows with tissue not in unique_tissues (ie. remove Normal tissue labeled rows)
subset_metadata = subset_metadata[subset_metadata['tissue'].isin(unique_tissues)]

# # group by Mouse and iterate over each group to write mouse specific data to files
# for mouse, df in subset_metadata.groupby('Mouse'):
#     mouse_outfile = outfile.replace(".csv", f"_{mouse}.csv")
#     df.to_csv(mouse_outfile, sep="\t", index=False)

# # write subset metadata to file
# subset_metadata.to_csv(outfile, sep="\t", index=False)



# repeat above but for ES clones to make sure a mouse with several tissues is not from a different source clone
es_clones = subset_metadata.groupby('Mouse')['ES_clone'].unique()



# get necessary lane IDs for the subset of mice for cassiopeia preprocessing of SRA data
lane_ids = subset_metadata['Lane'].unique().tolist()

# get necessary SRA IDs corresponding to the lane IDs
sra_metadata_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/targetsite_sraruninfo.csv"
outfile_sra_metadata = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/targetsite_sraruninfo_FILTERED.csv"

sra_metadata = pd.read_csv(sra_metadata_file, sep=",")
sra_metadata['Lane'] = sra_metadata['LibraryName'].str.split(" ").str[0]
subset_sra_metadata = sra_metadata[sra_metadata['Lane'].isin(lane_ids)]
sra_id_string = ""
for id in subset_sra_metadata['Run'].tolist():
    sra_id_string+=f"{id} "

# # write subset sra metadata to file
# subset_sra_metadata.to_csv(outfile_sra_metadata, sep="\t", index=False)


# output mouse specific allele files from the subset allele file
allele_filepath = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered.txt"
allele_outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/processed_metadata/KPTracer.alleleTable.FINAL.filtered.txt"
allele_df = pd.read_csv(allele_filepath, sep="\t")
allele_df = allele_df.drop(columns=['Unnamed: 0'])
allele_df["Mouse"] = allele_df['Tumor'].str.split("_").str[0]
# group by Mouse and iterate over each group to write mouse specific data to files
for mouse, df in allele_df.groupby('Mouse'):
    mouse_outfile = allele_outfile.replace(".txt", f"_{mouse}.txt")
    # df.to_csv(mouse_outfile, sep="\t", index=False)

    # get mouse specific cell tissue labels
    tissue_outfile = mouse_outfile.replace(".txt", f"_tissue_labels.txt")
    df['tissue'] = [''.join(filter(str.isalpha, name)) for name in df['Tumor'].str.split('_').str[2]]
    df = df[['cellBC', 'tissue']]
    df = df.drop_duplicates()
    df.to_csv(tissue_outfile, sep="\t", index=False)
    

# go from allele table to character matrix by combining intBC alleles
subset_allele_df = allele_df
intBCs_per_mouse = subset_allele_df.groupby('Mouse')['intBC'].unique().apply(len)

# write mouse specific character matrix to file
for mouse, df in subset_allele_df.groupby('Mouse'):

    # built in cassiopeia function to convert allele table to character matrix
    char_matrix = cas.pp.convert_alleletable_to_character_matrix(df)
    char_matrix_df = char_matrix[0]
    mut_dict = char_matrix[2]

    mouse_outfile = allele_outfile.replace(".txt", f"_{mouse}_char_matrix.txt")
    char_matrix_df.index.name = 'cellBC'
    char_matrix_df.to_csv(mouse_outfile, sep="\t", index=True)

    # write mutation dictionary to file
    mut_dict_outfile = mouse_outfile.replace(".txt", f"_mut_dict.txt")
    with open(mut_dict_outfile, "w") as f:
        f.write(f"site\tmut_id\tmut_str\n")
        for bc in mut_dict.keys():
            site = bc+1
            for key, value in mut_dict[bc].items():
                f.write(f"r{site}\t{key}\t{value}\n")
