
import sys
from Bio import Phylo

# user input
newick_file = sys.argv[1]

# format outfile
outfile = newick_file.replace(".nwk", ".nexus")

# read in newick tree
tree = Phylo.read(newick_file, "newick")

# write out nexus tree
Phylo.write(tree, outfile, "nexus")
