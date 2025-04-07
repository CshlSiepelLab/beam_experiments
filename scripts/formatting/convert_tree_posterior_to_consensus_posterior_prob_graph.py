#!/usr/bin/env python3

import sys
from beam_visualization import BeamResults


# user inputs
beam_trees = sys.argv[1]
beam_log = sys.argv[2]
primary_tissue = sys.argv[3]
outfile = sys.argv[4]
cores = int(sys.argv[5])
output_file_probability_plot = sys.argv[6]
thresholds = str(sys.argv[7]).split(",")
outprefix_thresholded = sys.argv[8]


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

# plot the probability consensus graph
results.plot_probability_graph(
    output_file=output_file_probability_plot,
)

# plot the probability consensus graph for several thresholds
results.plot_thresholded_graph(
    thresholds=thresholds,
    output_file_prefix=outprefix_thresholded,
)
