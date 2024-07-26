#!/usr/env/python3

import sys
import pandas as pd
import numpy as np
from ete3 import Tree

def get_site_category(label):
    site_category = ""
    if label == primary_tissue:
        site_category = "primary"
    else:
        site_category = "metastasis"
    return site_category

# tree=sys.argv[1]
# tissues=sys.argv[2]
# primary_tissue=sys.argv[3]
# outfile=sys.argv[4]

# for testing
treefile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/laml/2945/2945_laml_trees.nwk"
tissues="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/raw_data/2945/cell_tree_seed1082116693.labeling"
primary_tissue="P"
outdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/metient/2945"


# use tree to get edge list and branch lengths
tree = Tree(treefile, format=3)
tree.get_tree_root().name = '0'

edges = pd.DataFrame(columns = ['node1', 'node2'])
branch_lengths = {}
cluster_labels = []
cluster_index_label = {}

# set root info
cluster_index_label[tree.get_tree_root().name] = 0
branch_lengths[cluster_index_label[tree.get_tree_root().name]] = tree.get_tree_root().dist
cluster_labels.append(tree.get_tree_root().name)

# set all other node info by children traversal
i = 1
for node in tree.traverse():
    if node.is_leaf() == True:
        name = node.name
    else:
        node_name = node.name
        children = node.children
        for child in children:
            child_name = child.name
            cluster_labels.append(child_name)
            cluster_index_label[child_name] = i
            branch_lengths[i] = child.dist
            edges.loc[len(edges)] = [cluster_index_label[node_name], i]
            i = i+1

### output edges
edges.to_csv(f"{outdir}/tree.txt", sep=" ", index=False, header=False)


# read in tissues and format metadata tsv
tissues_df = pd.read_csv(tissues, sep=" ", names=["id", "tissue"])
tissues_dict = dict(zip(tissues_df["id"].astype(str), tissues_df["tissue"].astype(str)))
unique_tissues = set(tissues_dict.values())
tissue_to_int = {tissue: i for i, tissue in enumerate(unique_tissues)}

### output metadata tsv
with open(f"{outdir}/metadata.tsv", "w") as file:
    file.write("\t".join(["anatomical_site_index", "anatomical_site_label", "cluster_index", "cluster_label", "present", "site_category", "num_mutations"]))
    for cluster_label in cluster_labels:
        cluster_index = cluster_index_label[cluster_label]
        num_mutations = branch_lengths[cluster_index]
        for anatomical_site_label in unique_tissues:
            anatomical_site_index = tissue_to_int[anatomical_site_label]
            if anatomical_site_label != primary_tissue:
                site_category = "metastasis"
            else:
                site_category = "primary"
            if cluster_label in tissues_dict and tissues_dict[cluster_label] == anatomical_site_label:
                present = 1
            else:
                present = 0
            file.write(f"\n{anatomical_site_index}\t{anatomical_site_label}\t{cluster_index}\t{cluster_label}\t{present}\t{site_category}\t{num_mutations}")