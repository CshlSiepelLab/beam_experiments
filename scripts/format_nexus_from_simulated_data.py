#!/usr/bin/env python3

import sys
import ete3
import pandas as pd

def create_nexus_from_dict(data_dict):
    xml_content = ""

    for taxon_id, state in data_dict.items():
        xml_content += f"<sequence id=\"cell{taxon_id}\" spec=\"Sequence\" taxon=\"cell{taxon_id}\" value=\"{state[1:]}\"/>\n"

    return xml_content

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
xml_content = create_nexus_from_dict(leaf_tissues_dict)

# Output alignment seciton data for simulated samples
file_path = newick_filepath.split(".")[0] + "_xml_sample.xml"
with open(file_path, "w") as file:
    file.write(xml_content)

# Output tree in nexus file format
nexus_file_path = newick_filepath.split(".")[0] + ".nexus"
tree.write(format=9, outfile=nexus_file_path)