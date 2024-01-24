#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def replace_values_based_on_dictionary(df, column_dict):
    # Create a copy of the DataFrame to avoid modifying the original
    result_df = df.copy()
    # Iterate through the columns in the DataFrame
    for col in df.columns:
        # Check if the column is in the dictionary
        if col in column_dict:
            # Replace values of 1 with the corresponding unique integer from the dictionary
            result_df[col] = result_df[col].replace({1: column_dict[col]})
    return result_df

def partition_names_to_sites(names, sites):
    # Extract positions from names and convert to NumPy array
    positions = np.array([int(name.split('_')[1]) for name in names])
    # Find the closest site index for each position
    closest_sites = np.argmin(np.abs(positions[:, None] - np.array(sites)), axis=1)
    # Create the dictionary with names as keys and corresponding site index as values
    result_dict = {name: site_index for name, site_index in zip(names, closest_sites)}
    return result_dict

### User input parameters
evotracer_mut_matrix_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/MMUS1495_mutation_matrix.csv"
asv_stat_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/asv_stat_tissues.csv"
cutsites_str = "17,43,69,95,121,147,173,199,225,251"

og_matrix = pd.read_csv(evotracer_mut_matrix_path, index_col=0)
asv_tissues = pd.read_csv(asv_stat_path)
cutsites = cutsites_str.split(",")
cutsites = [int(site) for site in cutsites]

# Obtain all mutation names as np array
uniq_mutations = og_matrix.columns.values

# Keep original mutation names mapped to unique integer value for later post-processing
mut_dict = {name : (i+1) for i, name in enumerate(uniq_mutations)}

# Replace binary presence with unique integer for the mutations
int_matrix = replace_values_based_on_dictionary(og_matrix, mut_dict)

# Partition mutation names to site with sites 0 indexed
mut_sites_dict = partition_names_to_sites(uniq_mutations, cutsites)

# Rename column names based on site and collapse same sites with comma sep between more than one value for the same site to then give error message if this occurs
site_matrix = int_matrix.rename(columns=mut_sites_dict, inplace=False)
grouped_site_matrix = site_matrix.groupby(site_matrix.columns, axis=1)
for group in grouped_site_matrix:
    site_num = group[0]
    num_mutations = len(group[1].columns.values)
    print(f"Site {site_num} has {num_mutations} unique mutation(s).")
collapsed_matrix = grouped_site_matrix.agg(lambda x: ','.join(str(val) for val in x[x != 0]))
