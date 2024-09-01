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
        print(name)
        tissue = tissue_df.loc[tissue_df['node'] == name, 'tissue'].values[0]
        traits += f"{name}={tissue}"
        if i < len(taxa_names) - 1:
            traits += ",\n"
        else:
            traits += "\n"
    return traits

newick_file = sys.argv[1]
tissue_file = sys.argv[2]

# newick_file = "gundem_a10/A10_unlabeled_tree.nwk"
# tissue_file = "gundem_a10/A10_tissues.tsv"

# newick_file = "results/unsymmetrical_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_12_24/machina_m8_sim_data/seed0/T_seed0_unlabeled_true_tree.nwk"
# tissue_file = "results/unsymmetrical_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_12_24/machina_m8_sim_data/seed0/T_seed0_tissues.tsv"

try:
    tree = Tree(newick_file, format=5)
except:
    tree = Tree(newick_file, format=8)
tissue_df = pd.read_csv(tissue_file, sep='\t')
tissue_df = tissue_df.loc[:, ['node', 'tissue']]

# Rename newick string names with cell prefix for consistency if any name is all int values
taxa_names = []
add_cell = False
for leaf in tree.iter_leaves():
    current_name = leaf.name
    try:
        new_name = "cell" + int(current_name)
        add_cell = True
    except ValueError:
        new_name = current_name
    leaf.name = new_name
    taxa_names.append(new_name)

if add_cell:
    tissue_df['node'] = 'cell' + tissue_df['node'].astype(str)

# Replace semiccolons from machina sims; should not affect my own sim data
tissue_df['node'] = tissue_df['node'].str.replace(";", "_")

# Output relabeled newick string
newick_outfile = newick_file.split(".")[0] + "_newick_formatted_for_xml.txt"
newick = tree.write(format=5, format_root_node=False)
# removes outer parentheses to set unedited as root length when fake empty root exists as artifact from ete3 tree building
if len(tree.get_tree_root().children) == 1:
    newick = newick[1:-2] + "\n"
else:
    newick = newick.replace(";", "") + "\n"
    
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
