
import sys

from graphposterior.format_inputs import format_beam_inputs


indel_matrix_file = sys.argv[1]
outdir = sys.argv[2]

format_beam_inputs(indel_matrix_file, outdir)

