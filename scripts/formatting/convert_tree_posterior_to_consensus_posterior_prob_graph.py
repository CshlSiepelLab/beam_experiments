#!/usr/bin/env python3

import sys
from beam_visualization import BeamResults


# user inputs
beam_trees = sys.argv[1]
beam_log = sys.argv[2]
primary_tissue = sys.argv[3]
outfile = sys.argv[4]
cores = int(sys.argv[5])

# load in data
results = BeamResults(
            beam_trees, 
            beam_log, 
            primary_tissue=primary_tissue
            )

# get consensus graph and write to file
results.get_consensus_graph(
    cores=cores,
    output_file=outfile
)

