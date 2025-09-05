
import sys
import os
import pandas as pd

from beam_sup.matrix_utils import expand_clones_multiple_tissues

# inputs
matrix_file = sys.argv[1]
tissues_file = sys.argv[2]
output_dir = sys.argv[3]

# read in matrix and tissues files to df where the files are tsv files
matrix_df = pd.read_csv(matrix_file, sep="\t", index_col=0)
tissues_df = pd.read_csv(tissues_file, sep="\t")

matrix_df, tissues_df = expand_clones_multiple_tissues(matrix_df, tissues_df)

# write outputs
matrix_df.to_csv(os.path.join(output_dir, "expanded_clones_matrix.tsv"), sep="\t")
tissues_df.to_csv(
    os.path.join(output_dir, "expanded_clones_tissues.tsv"), sep="\t", index=False
)
