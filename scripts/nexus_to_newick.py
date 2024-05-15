#!/usr/bin/env python3

import sys
from Bio import Phylo

# user input 
nexus_file = sys.argv[1]

# # testing
# nexus_file="results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/mS/20692/joint_inference_beast_tissues.tree"

# format outfile
outfile=nexus_file.replace(".tree",".nwk")

# read in tree and print as newick
tree = Phylo.read(nexus_file, 'nexus')
Phylo.write(tree, outfile, 'newick')
