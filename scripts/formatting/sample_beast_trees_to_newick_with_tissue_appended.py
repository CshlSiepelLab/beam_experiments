
import sys
import os

from beam_sup.posterior_processing import sample_beast_trees_and_append_tissue_to_names

nexus_file = sys.argv[1]

outdir = os.path.dirname(nexus_file)

sample_beast_trees_and_append_tissue_to_names(nexus_file, outdir, 10)