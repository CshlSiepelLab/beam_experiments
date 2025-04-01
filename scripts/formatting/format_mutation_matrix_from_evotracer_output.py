#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import pdb


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
    positions = np.array([int(name.split("_")[1]) for name in names])
    # Find the closest site index for each position
    closest_sites = np.argmin(np.abs(positions[:, None] - np.array(sites)), axis=1)
    # Create the dictionary with names as keys and corresponding site index as values
    result_dict = {name: site_index for name, site_index in zip(names, closest_sites)}
    return result_dict


def aggregate_site_matrix(grouped_site_matrix):
    group_names = [group[0] for group in grouped_site_matrix]
    asv_names = grouped_site_matrix.get_group(0).index.values
    agg_df = pd.DataFrame(index=asv_names, columns=group_names)
    for group in grouped_site_matrix:
        group_name = group[0]
        for index, row in group[1].iterrows():
            row_name = row.name
            row_concat = ",".join([str(val) for val in row if val != 0])
            if row_concat == "":
                row_concat = "0"
            agg_df.loc[row_name, group_name] = row_concat
    return agg_df


def replace_multiple_mutations(agg_df, mut_dict):
    new_dict = mut_dict.copy()
    joint_mut_df = agg_df.copy()
    all_values = list(joint_mut_df.values.flatten())
    current_key = int(list(new_dict.values())[-1]) + 1
    column_names = list(joint_mut_df.columns)
    for index, row in joint_mut_df.iterrows():
        for col_name in column_names:
            value = joint_mut_df.loc[index, col_name]
            count = all_values.count(value)
            if "," in str(value) and value not in new_dict and count > 1:
                joint_mut_df.loc[index, col_name] = current_key
                new_dict[value] = current_key
                current_key += 1
            elif "," in str(value) and count > 1:
                joint_mut_df.loc[index, col_name] = new_dict[value]
            elif "," in str(value):
                joint_mut_df.loc[index, col_name] = "0"
    return joint_mut_df, new_dict


### User input parameters
evotracer_mut_matrix_path = sys.argv[1]
# evotracer_mut_matrix_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/raw_data/MMUS1495_mutation_matrix.csv"
# asv_stat_path = "/Users/staklins/projects/crispr-barcode-cancer-metastasis/bayesian_phylogenetic_metastasis/examples/real_data/mmus1495/raw_data/asv_stat_tissues.csv"
cutsites_str = "17,43,69,95,121,147,173,199,225,251"

og_matrix = pd.read_csv(evotracer_mut_matrix_path, index_col=0)
# asv_tissues = pd.read_csv(asv_stat_path)
cutsites = cutsites_str.split(",")
cutsites = [int(site) for site in cutsites]

# Count non-zero values in each column
non_zero_counts = og_matrix.astype(bool).sum(axis=0)

# Drop columns with 0 or 1 count of non-zero values to remove non-present mutations and singleton mutations
columns_to_drop = non_zero_counts[non_zero_counts <= 1].index
og_matrix = og_matrix.drop(columns=columns_to_drop)

# Obtain all mutation names as np array
uniq_mutations = og_matrix.columns.values

# Keep original mutation names mapped to unique integer value for later post-processing
mut_dict = {name: (i + 1) for i, name in enumerate(uniq_mutations)}

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
    count_list = group[1].apply(lambda row: (row != 0).sum(), axis=1).values
    count_above_one = sum(1 for value in count_list if value > 1)
    print(
        f"Site {site_num} has {count_above_one} occurences of more than 1 mutation in the same ASV."
    )

collapsed_matrix = aggregate_site_matrix(grouped_site_matrix)
joint_mut_matrix, joint_mut_dict = replace_multiple_mutations(
    collapsed_matrix, mut_dict
)

output_path = evotracer_mut_matrix_path.split(".")[0] + "_reformatted_tidetree.tsv"
joint_mut_matrix.to_csv(output_path, sep="\t")

output_dict_path = (
    evotracer_mut_matrix_path.split(".")[0] + "_reformatted_tidetree_dictionary.tsv"
)
with open(output_dict_path, "w") as file:
    file.write(f"original_mutation_name\ttidetree_mutation_int\n")
    for key, value in joint_mut_dict.items():
        file.write(f"{key}\t{value}\n")
