#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def format_sequences_from_matrix(mut_df):
    seqs = np.empty(len(mut_df), dtype=object)
    mut_dict = {0:0}
    i = 1
    
    for index, values in mut_df.iterrows():
        row_index = values.values[0]
        raw_values = values.values[1:].astype(int)
        unique_values = np.unique(raw_values)
        unique_values = unique_values[unique_values != 0]
        for value in unique_values:
            if value not in mut_dict.keys():
                mut_dict[value] = i
                i = i+1
        # Replacing unique integers indicating mutations in the matrix with consecutive unique integers for tidetree compatibility
        new_values = np.array([mut_dict[key] for key in raw_values])
        row_concat = ','.join(map(str, new_values))
        new_string = f"<sequence id='cell{row_index}' spec='Sequence' taxon='{index+1}' value='{row_concat},'/>\n"
        seqs[index] = new_string
    return seqs, mut_dict

### User input parameters
mutation_matrix_filepath=sys.argv[1]
#mutation_matrix_filepath="examples/simulated_data/sim_results_ten_samples/ten_samples_indel_character_matrix.tsv"

matrix_df = pd.read_csv(mutation_matrix_filepath, sep="\t")
sequences, tidetree_dict = format_sequences_from_matrix(matrix_df)

# Output xml formatted strings
output_path = mutation_matrix_filepath.split(".")[0] + "_tidetree_sequences.xml"
with open(output_path, "w") as file:
    for seq in sequences:
        file.write(seq)
        
# Output mutation mapping of integer representations
dict_df = pd.DataFrame(list(tidetree_dict.items()), columns=['raw_mutation_number', 'tidetree_replacement'])
output_path = mutation_matrix_filepath.split(".")[0] + "_tidetree_mutation_dict.tsv"
dict_df.to_csv(output_path, sep = '\t', index=False)
