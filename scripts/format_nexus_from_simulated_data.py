#!/usr/bin/env python3

import sys
import ete3
import pandas as pd

def create_nexus_from_dict(data_dict, site_number):
    nexus_content = "#NEXUS\nBEGIN DATA;\n"
    
    # Determine the number of taxa
    taxa_count = len(data_dict)
    
    # Write the dimensions and format
    nexus_content += f"DIMENSIONS NTAX={taxa_count} NCHAR={site_number};\n"
    nexus_content += "FORMAT DATATYPE=STANDARD MISSING=? GAP=-;\n"
    
    # Write the matrix
    nexus_content += "MATRIX\n"
    
    for taxon_id, state in data_dict.items():
        # Use taxon_id as the taxon label
        ### TEMP FIX TO ONLY LABEL NUMERICALLY
        nexus_content += f"{taxon_id} {state[1]}\n"
    
    nexus_content += ";\n"
    nexus_content += "END;\n"
    
    return nexus_content

# newick_filepath = str(sys.argv[1])
# tissue_data_filepath = str(sys.argv[2])
newick_filepath = "simulated_data/sim_results_test_sim/test_sim_true.nwk"
tissue_data_filepath = "simulated_data/sim_results_test_sim/test_sim_tissues.tsv"

tree = ete3.Tree(newick_filepath)
tissues_df = pd.read_csv(tissue_data_filepath, sep="\t")

leaf_names = [int(leaf.name) for leaf in tree.iter_leaves()]

leaf_tissues_df = tissues_df[tissues_df['node'].isin(leaf_names)]

leaf_tissues_dict = dict(zip(leaf_tissues_df['node'], leaf_tissues_df['tissue']))

# Make nexus file for a single site of tissue label aligned across samples
nexus_content = create_nexus_from_dict(leaf_tissues_dict, 1)

file_path = newick_filepath.split(".")[0] + ".nex"
with open(file_path, "w") as nexus_file:
    nexus_file.write(nexus_content)
