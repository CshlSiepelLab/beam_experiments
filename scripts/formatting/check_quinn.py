import os
from glob import glob
import pandas as pd

# Basic script to make sure my re-doing of the Quinn et al. 2021 preprocessing matches the deposited data in terms of number of cells and sites per CP

deposited_data_dir = "/grid/siepel/home/staklins/projects/crispr_barcode/data/quinn_2021_real_data/GSE161363/character_matrices"
my_processed_data_dir = "/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/quinn_2021_lung_cancer_data/successive_raw_data/5k"

deposited_char_matrix_files = glob(f"{deposited_data_dir}/*_matrix.alleleThresh.txt")
my_char_matrix_files = glob(f"{my_processed_data_dir}/*_successive_character_matrix.tsv")

num_matching_files = 0
num_mismatched_files = 0
for file in my_char_matrix_files:
    cp_num = os.path.basename(file).split("_")[0]
    my_char_matrix = pd.read_csv(file, sep="\t", index_col=0)
    matching_deposited_files = [f for f in deposited_char_matrix_files if f"lg{cp_num }_" in os.path.basename(f)]
    assert len(matching_deposited_files) == 1, f"Expected one matching deposited file for clone {cp_num}, found {len(matching_deposited_files)}"
    deposited_char_matrix = pd.read_csv(matching_deposited_files[0], sep="\t", index_col=0)
    
    if my_char_matrix.shape != deposited_char_matrix.shape:
        print(f"Shape mismatch for clone {cp_num}: my shape {my_char_matrix.shape}, deposited shape {deposited_char_matrix.shape}")
        num_mismatched_files += 1
    else:
        num_matching_files += 1

print(f"Number of matching files: {num_matching_files}")
print(f"Number of mismatched files: {num_mismatched_files}")

