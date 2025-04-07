#!/usr/bin/env python3

import sys
from beam_visualization import BeamResults


beam_trees = sys.argv[1]
beam_log = sys.argv[2]
primary_tissue = sys.argv[3]
output_file_matrix = sys.argv[4]
output_file_information = sys.argv[5]
threads = int(sys.argv[6])


# Load in data to BeamResults object
results = BeamResults(beam_trees, beam_log, primary_tissue)

# Get migration count matrix and mutual information
results.compute_posterior_mutual_info(
    output_file_matrix = output_file_matrix,
    output_file_information = output_file_information,
    threads = threads
)

