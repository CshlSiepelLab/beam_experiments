
import sys
import ete3
import pandas as pd
from collections import Counter
import random

from beam_sup.migration_inference import label_tissues_fitch_parsimony, label_tissues_consensus, label_tissues_random


tree_file = sys.argv[1]
leaf_tissues_tsv = sys.argv[2]
outdir = sys.argv[3]

primary_tissue = "P"

tree = ete3.Tree(tree_file, format=8)

# Force the tree to be binary (ie. resolve polytomies)
tree.resolve_polytomy()

tissues_df = pd.read_csv(
    leaf_tissues_tsv,
    sep=r"\s+",
    header=None,
    names=["cell", "tissue"],
    dtype={"cell": str, "tissue": str},
)

# Get migration inference results
random_tree = label_tissues_random(tree, 
                                   tissues_df, 
                                   primary_tissue)
consensus_tree = label_tissues_consensus(tree, 
                                         tissues_df, 
                                         primary_tissue)
parsimony_tree, all_solutions = label_tissues_fitch_parsimony(tree, 
                                                              tissues_df, 
                                                              primary_tissue, 
                                                              1)

# Output trees to newick files
random_output = outdir + "/random_tissues.nwk"
random_tree.write(outfile=random_output, format=8, format_root_node=True)

consensus_output = outdir + "/consensus_tissues.nwk"
consensus_tree.write(outfile=consensus_output, format=8, format_root_node=True)

parsimony_output = outdir + "/parsimony_tissues.nwk"
parsimony_tree.write(outfile=parsimony_output, format=8, format_root_node=True)
