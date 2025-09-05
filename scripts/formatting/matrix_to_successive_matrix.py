
import sys
import os
import pandas as pd

from beam_sup.matrix_utils import convert_matrix_to_row_successive_matrix

char_matrix_file = sys.argv[1]
mut_dict_file = sys.argv[2]
outdir = sys.argv[3]

char_matrix_df = pd.read_csv(char_matrix_file, sep="\t", index_col=0)
mut_dict_df = pd.read_csv(mut_dict_file, sep="\t", index_col=0, header=None)
mut_dict = {row[1][0]: row[0] for row in mut_dict_df.iterrows()}

successive_char_matrix, successive_mut_dict, _ = convert_matrix_to_row_successive_matrix(char_matrix_df, mut_dict)

outfile_char_matrix = f"{outdir}/successive_char_matrix.csv"
outfile_mut_dict = f"{outdir}/successive_mut_dict.csv"

successive_char_matrix.to_csv(outfile_char_matrix, sep="\t")
with open(outfile_mut_dict, "w") as f:
    f.write("mut_id\tmut_str\n")
    for mut_str, mut_id in successive_mut_dict.items():
        f.write(f"{mut_id}\t{mut_str}\n")
