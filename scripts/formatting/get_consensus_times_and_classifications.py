
import sys

from beam.posterior_processing import get_consensus_times_and_classifications


file_path = sys.argv[1]
origin_tissue = sys.argv[2]
consensus_threshold = float(sys.argv[3])
outprefix = sys.argv[4]

get_consensus_times_and_classifications(file_path, origin_tissue, consensus_threshold, outprefix)