#!/usr/bin/env python3

import sys
import pandas as pd

# user input
# mutation_matrix_file = sys.arv[1]
# tissue_map_file = sys.argv[2]

# testing
mutation_matrix_file = "results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/1983/1983_indel_character_matrix.tsv"
tissue_map_file = "results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/1983/cell_tree_seed1983.labeling"

# make outfile
outfile = mutation_matrix_file.replace(".tsv", "_with_tissues.tsv")

# read in files to df
mutation_matrix = pd.read_csv(mutation_matrix_file, sep="\t", header=0, index_col=0)
tissue_map = pd.read_csv(tissue_map_file, sep=" ", names=["tissue"], index_col=0)

# append tissues to mutation matrix
tissue_mut_matrix = pd.concat([mutation_matrix, tissue_map["tissue"]], axis=1)

# write new mutation matrix with tissue to tsv
tissue_mut_matrix.to_csv(outfile, sep="\t")
