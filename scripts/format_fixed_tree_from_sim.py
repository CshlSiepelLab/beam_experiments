#!/usr/bin/env python3

import sys
import ete3
import pandas as pd

def format_tree_file(data_dict):
    xml_content = f"#NEXUS\n\nBegin taxa;\n\tDimensions ntax={len(data_dict)};\n\t\tTaxlabels\n"
    translate_content = "Begin trees;\n\tTranslate\n"
    translate_dict = {}

    for index, (taxon_id, state) in enumerate(data_dict.items()):
        # xml_content += f"<sequence id=\"{taxon_id}\" spec=\"Sequence\" taxon=\"{taxon_id}\" value=\"{state[1:]}\"/>\n"
        xml_content += "\t\t\tcell" + str(taxon_id) + "\n"
        translate_content += f"\t\t\t{index + 1} cell{taxon_id},\n"
        translate_dict[str(taxon_id)] = index + 1   # Need to index from 1 for BEAUti
    xml_content = xml_content + "\t\t\t;\nEnd;\n" + translate_content[:-2] + "\n\t\t\t;\n"

    return xml_content, translate_dict

newick_filepath = str(sys.argv[1])
tissue_data_filepath = str(sys.argv[2])

tree = ete3.Tree(newick_filepath)
tissues_df = pd.read_csv(tissue_data_filepath, sep="\t")

leaf_names = [int(leaf.name) for leaf in tree.iter_leaves()]
leaf_tissues_df = tissues_df[tissues_df.loc[:,'node'].isin(leaf_names)]
leaf_tissues_dict = dict(zip(leaf_tissues_df.loc[:,'node'], leaf_tissues_df.loc[:, 'tissue']))

# Make nexus file for a single site of tissue label aligned across samples
xml_content, translate_dict = format_tree_file(leaf_tissues_dict)

# Replace newick keys with translate value indices from 1 up to total number of taxa (requirement of BEAUti to load in trees)
for leaf in tree.iter_leaves():
    leaf.name = translate_dict[leaf.name]

# Output tree in nexus file format
### Hack to prevent beast error for extra node at root; Need to solve this another way eventually
# Use below to print both leaf names and all branch lengths for the tree
nexus_file_path = newick_filepath.split(".")[0] + ".tree"
newick = tree.write(format=5)
indexes = [index for index,char in enumerate(newick) if char == ":"]
index = indexes[-1]
newick_fixed = newick[1:index] + ";"

# Use below to print onlt the leaf names of tree
# newick = tree.write(format=9)
# indexes = [index for index,char in enumerate(newick) if char == ")"]
# index = indexes[-1]
# newick_fixed = newick[1:index] + ";"


xml_content = xml_content + "tree TREE1 = " + newick_fixed + "\nEnd;"
with open(nexus_file_path, "w") as file:
    file.write(xml_content)

# Output tissues data in tsv
tissues_path = newick_filepath.split(".")[0] + ".dat"
tissues_key = leaf_tissues_df[['node','tissue']]
tissues_key.loc[:, 'node'] = 'cell' + tissues_key.loc[:, 'node'].astype(str)
tissues_key.to_csv(tissues_path, sep="\t", index=False, header=False)