#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def format_sequences_from_matrix(mut_df):
    seqs = np.empty(len(mut_df), dtype=object)
    for index, values in mut_df.iterrows():
        row_index = values.values[0]
        row_concat = ','.join(map(str, values.values[1:]))
        new_string = f"<sequence id='cell{row_index}' spec='Sequence' taxon='{index}' value='{row_concat},'/>\n"
        seqs[index] = new_string
    return seqs

### User input parameters
mutation_matrix_filepath=sys.argv[1]
#mutation_matrix_filepath="examples/simulated_data/sim_results_ten_samples/ten_samples_indel_character_matrix.tsv"

matrix_df = pd.read_csv(mutation_matrix_filepath, sep="\t")
sequences = format_sequences_from_matrix(matrix_df)

output_path = mutation_matrix_filepath.split(".")[0] + ".xml"
with open(output_path, "w") as file:
    for seq in sequences:
        file.write(seq)
