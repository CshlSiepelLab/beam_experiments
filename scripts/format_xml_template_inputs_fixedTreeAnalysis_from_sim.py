#!/usr/bin/env python3

# This script takes in a newick file and a tsv file with tissue label mapping for the names in the newick file. It outputs two files, one with the sequence input and the other with the trait input for the xml file formatting to run FixedTreeAnalysis.

import sys
import pandas as pd
from ete3 import Tree

def format_sequences_string(name_list):
    sequences = ""
    i = 0
    for name in name_list:
        new_string = f"<sequence id='Sequence.{i}' spec='Sequence' taxon='{name}' totalcount='4' value='?'/>\n"
        sequences += new_string
        i = i + 1
    return sequences

def format_taxa_string(name_list):
    taxa = ""
    for name in name_list:
        taxa += f"<taxon id='{name}' spec='Taxon'/>\n"
    return taxa

def format_traitset(taxa_names,tissue_df):
    traits = ""
    for i, name in enumerate(taxa_names):
        tissue = tissue_df.loc[tissue_df['node'] == name, 'tissue'].values[0]
        traits += f"{name}={tissue}"
        if i < len(taxa_names) - 1:
            traits += ",\n"
        else:
            traits += "\n"
    return traits

newick_file = sys.argv[1]
tissue_file = sys.argv[2]

# newick_file = "compare_beast_machina_fixedtree_2_1_24/sim_results_sim1/sim1_true.nwk"
# tissue_file = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/compare_beast_machina_fixedtree_2_1_24/sim_results_sim1/sim1_tissues.tsv"

tree = Tree(newick_file, format=5)
tissue_df = pd.read_csv(tissue_file, sep='\t')
tissue_df = tissue_df.loc[:, ['node', 'tissue']]

# Rename newick string names with cell prefix for consistency
taxa_names = []
for leaf in tree.iter_leaves():
    current_name = leaf.name
    new_name = "cell" + current_name
    leaf.name = new_name
    taxa_names.append(new_name)

tissue_df['node'] = 'cell' + tissue_df['node'].astype(str)

# Output relabeled newick string
newick_outfile = newick_file.split(".")[0] + "_newick_formatted_for_xml.txt"
newick = tree.write(format=5, format_root_node=False)
# removes outer parentheses to set unedited as root length
newick = newick[1:-2] + "\n"
with open(newick_outfile, "w") as file:
    file.write(newick)
    
# Output sequence data
sequences_outfile = newick_file.split(".")[0] + "_sequences_formatted_for_xml.txt"
sequences = format_sequences_string(taxa_names)
with open(sequences_outfile, "w") as file:
    file.write(sequences)

# Output taxon set
taxon_outfile = newick_file.split(".")[0] + "_taxonset_formatted_for_xml.txt"
taxa = format_taxa_string(taxa_names)
with open(taxon_outfile, "w") as file:
    file.write(taxa)

# Output trait set
trait_outfile = newick_file.split(".")[0] + "_traitset_formatted_for_xml.txt"
traits = format_traitset(taxa_names,tissue_df)
with open(trait_outfile, "w") as file:
    file.write(traits)
