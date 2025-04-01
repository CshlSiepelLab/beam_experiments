#!/usr/bin/env python3

# This script takes in a newick tree file and a csv file with leaf names matching the newick in the first column and the tissue labels in the second column. Multiple rows with the same leaf name in the input csv file are allowed and used to make duplicate sister leaves so that there is a 1 to 1 mapping of leaves to sampled tissues.

import sys
import ete3
import pandas as pd


def format_tree_file(data_dict):
    xml_content = (
        f"#NEXUS\n\nBegin taxa;\n\tDimensions ntax={len(data_dict)+1};\n\t\tTaxlabels\n"
    )
    translate_content = "Begin trees;\n\tTranslate\n"
    translate_dict = {}
    for index, (taxon_id, state) in enumerate(data_dict.items()):
        xml_content += "\t\t\t" + str(taxon_id) + "\n"
        max_val = index + 1
        translate_content += f"\t\t\t{max_val} {taxon_id},\n"
        translate_dict[str(taxon_id)] = max_val  # Need to index from 1 for BEAUti
    translate_dict["Root"] = max_val + 1
    xml_content += "\t\t\t" + "Root" + "\n"
    translate_content += f"\t\t\t{max_val + 1} Root,\n"
    xml_content = (
        xml_content + "\t\t\t;\nEnd;\n" + translate_content[:-2] + "\n\t\t\t;\n"
    )
    return xml_content, translate_dict


newick_filepath = str(sys.argv[1])
tissue_data_filepath = str(sys.argv[2])
root_trait = "PRL"

tree = ete3.Tree(newick_filepath, format=9)
tissues_df = pd.read_csv(tissue_data_filepath, sep=",")

leaf_names = [str(leaf.name) for leaf in tree.iter_leaves()]
leaf_tissues_df = tissues_df[tissues_df.loc[:, "node"].isin(leaf_names)]
leaf_tissues_df = leaf_tissues_df.sort_values("node")

counts = leaf_tissues_df["node"].value_counts()

for leaf in tree.iter_leaves():
    name = leaf.name
    total = counts[name]
    if total > 1:
        i = 1
        while i != total:
            clone = leaf.copy()
            clone.name = clone.name + f"-{i}"
            leaf.add_sister(clone)
            i = i + 1

leaf_tissues_dict = {}
for index, (key, value) in leaf_tissues_df.iterrows():
    if key in leaf_tissues_dict:
        new_key = f"{key}-{j}"
        leaf_tissues_dict[new_key] = value
        j = j + 1
    else:
        leaf_tissues_dict[key] = value
        j = 1

# Make nexus file for a single site of tissue label aligned across samples
xml_content, translate_dict = format_tree_file(leaf_tissues_dict)

for node in tree.traverse():
    if node.is_leaf() is False:
        node.name = ""

# Replace newick keys with translate value indices from 1 up to total number of taxa (requirement of BEAUti to load in trees)
for leaf in tree.iter_leaves():
    leaf.name = translate_dict[leaf.name]

# set root to new index for fixed root trait analysis
root = tree.get_tree_root()
root.name = translate_dict["Root"]
root.dist = 0

leaf_tissues_dict["Root"] = root_trait

# Output tree in nexus file format
nexus_file_path = newick_filepath.split(".")[0] + ".tree"
newick = tree.write(format=5, format_root_node=True)

xml_content = (
    xml_content
    + "tree TREE1 = ("
    + newick[:-1]
    + f",{translate_dict['Root']}:0);\nEnd;"
)
with open(nexus_file_path, "w") as file:
    file.write(xml_content)

# Output tissues data in tsv
tissues_path = newick_filepath.split(".")[0] + ".dat"
with open(tissues_path, "w") as file:
    for key, value in leaf_tissues_dict.items():
        file.write(f"{key}\t{value}\n")
