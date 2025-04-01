#!/usr/bin/env python3

import sys
from Bio import Phylo

# # user input
# newick_file = sys.argv[1]

# testing
newick_file = "seed0.nwk"

# format outfile
outfile = newick_file.replace(".nwk", ".nexus")

# read in newick tree
tree = Phylo.read(newick_file, "newick")

# write out nexus tree
Phylo.write(tree, outfile, "nexus")
