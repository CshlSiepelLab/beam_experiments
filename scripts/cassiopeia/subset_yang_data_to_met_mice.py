#!/usr/bin/env python3

import pandas as pd

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

# repeat above but for ES clones to make sure a mouse with several tissues is not from a different source clone
es_clones = subset_metadata.groupby('Mouse')['ES_clone'].unique()


# group by Mouse and iterate over each group to write mouse specific data to files
for mouse, df in subset_metadata.groupby('Mouse'):
    mouse_outfile = outfile.replace(".csv", f"_{mouse}.csv")
    df.to_csv(mouse_outfile, sep="\t", index=False)

# # write subset metadata to file
# subset_metadata.to_csv(outfile, sep="\t", index=False)
