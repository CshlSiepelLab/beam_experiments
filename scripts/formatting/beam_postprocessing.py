
import sys
from graphposterior import PosteriorResults

import time


# user inputs
beam_trees = sys.argv[1]
beam_log = sys.argv[2]
primary_tissue = sys.argv[3]
outfile = sys.argv[4]
cores = int(sys.argv[5])
output_file_probability_plot = sys.argv[6]
thresholds = str(sys.argv[7]).split(",")
outprefix_thresholded = sys.argv[8]
outprefix_samples = sys.argv[9]
num_samples = int(sys.argv[10])
origin_time = int(sys.argv[11])
output_file_matrix = sys.argv[12]
output_file_information = sys.argv[13]
outprefix_metastasis_timing = sys.argv[14]
consensus_timing_threshold = float(sys.argv[15])

start_time = time.time()

# load in data
results = PosteriorResults(trees_file=beam_trees, log_file=beam_log, primary_tissue=primary_tissue, total_time=origin_time, cores=cores)
print("Loaded in data. Time elapsed: %.2f seconds" % (time.time() - start_time))

# get consensus graph and write to file
results.get_consensus_graph(output_file=outfile)
print("Computed consensus graph. Time elapsed: %.2f seconds" % (time.time() - start_time))

# plot the probability consensus graph
results.plot_probability_graph(output_file=output_file_probability_plot)
print("Plotted probability graph. Time elapsed: %.2f seconds" % (time.time() - start_time))

# plot the probability consensus graph for several thresholds
results.plot_thresholded_graph(threshold=thresholds, output_file_prefix=outprefix_thresholded)
print("Plotted thresholded graphs. Time elapsed: %.2f seconds" % (time.time() - start_time))

# sample trees from the posterior and plot
results.sample_and_plot_trees(n=num_samples, output_prefix=outprefix_samples)
print("Sampled and plotted trees. Time elapsed: %.2f seconds" % (time.time() - start_time))

# Get migration count matrix and mutual information
results.compute_posterior_mutual_info(output_file_matrix = output_file_matrix, output_file_information = output_file_information)
print("Computed migration matrix and mutual information. Time elapsed: %.2f seconds" % (time.time() - start_time))

# Get metastasis timing plot
results.get_metastasis_times(output_prefix = outprefix_metastasis_timing, min_prob_threshold= consensus_timing_threshold)
print("Computed metastasis timing. Time elapsed: %.2f seconds" % (time.time() - start_time))
