
import sys

from beam_sup.format_inputs import get_fixed_edit_rates_for_beam_sim_matrix

original_rates_file = sys.argv[1]
reordering_dict_file = sys.argv[2]
outfile = sys.argv[3]

get_fixed_edit_rates_for_beam_sim_matrix(
    original_rates_file, reordering_dict_file, outfile
)
