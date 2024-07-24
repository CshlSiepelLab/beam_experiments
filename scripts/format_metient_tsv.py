#!/usr/env/python3

import sys
import pandas as pd
import numpy as np

def get_site_category(label):
    site_category = ""
    if label == primary_tissue:
        site_category = "primary"
    else:
        site_category = "metastasis"
    return site_category

indel_matrix=sys.argv[1]
tissues=sys.argv[2]
primary_tissue=sys.argv[3]
outfile=sys.argv[4]

# # for testing
# indel_matrix="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/raw_data/2945/2945_indel_character_matrix.tsv"
# tissues="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/raw_data/2945/cell_tree_seed1082116693.labeling"
# primary_tissue="P"
# outfile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/metient/2945/metadata.tsv"

# read in tsv files
indel_df = pd.read_csv(indel_matrix, sep="\t", index_col=0)
indel_df.index = indel_df.index.astype(str)
tissues_df = pd.read_csv(tissues, sep=" ", names=["id", "tissue"])

# make df for output tsv
columns = ["anatomical_site_index", "anatomical_site_label", "cluster_index", "cluster_label", "present", "site_category", "num_mutations"]
df = pd.DataFrame(columns=columns)

unique_tissues = tissues_df["tissue"].unique().tolist()
tissue_to_int = {tissue: i for i, tissue in enumerate(unique_tissues)}

unique_ids = tissues_df['id']
id_to_int = {id: i for i, id in enumerate(unique_ids)}

entries = []

for id in tissues_df["id"]:
    anatomical_site_label = tissues_df[tissues_df['id'] == id]['tissue'].iloc[0]
    anatomical_site_index = tissue_to_int[anatomical_site_label]
    cluster_label = id
    cluster_index = id_to_int[id]
    present = 1
    site_category = get_site_category(anatomical_site_label)
    num_mutations = np.sum(indel_df.loc[str(id)].values > 0)

    row = {"anatomical_site_index": anatomical_site_index, 
           "anatomical_site_label": anatomical_site_label, 
           "cluster_index": cluster_index, 
           "cluster_label": cluster_label, 
           "present": present, 
           "site_category": site_category, 
           "num_mutations": num_mutations}
    entries.append(row)

    # Add rows for cell absent in other sites since all id's are single cells and metient needs each cluster to have a row for each site
    missing_sites = [site for site in unique_tissues if site != anatomical_site_label]
    for site in missing_sites:
        anatomical_site_label = site
        anatomical_site_index = tissue_to_int[anatomical_site_label]
        present = 0
        site_category = get_site_category(site)
    
        row = {"anatomical_site_index": anatomical_site_index, 
            "anatomical_site_label": anatomical_site_label, 
            "cluster_index": cluster_index, 
            "cluster_label": cluster_label, 
            "present": present, 
            "site_category": site_category, 
            "num_mutations": num_mutations}
        entries.append(row)
        

df = pd.DataFrame(entries)

df = df.sort_values(by=["cluster_index", "anatomical_site_index"])

df.to_csv(outfile, sep="\t", index=False)

