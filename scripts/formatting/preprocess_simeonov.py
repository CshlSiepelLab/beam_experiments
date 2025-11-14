
from glob import glob
import os
import pandas as pd

from beam_sup.matrix_utils import count_informative_characters, expand_clones_with_multiple_tissues


def convert_simeonov_barcode_file_to_indel_matrix(input_file, output_file):
    """
    Converts a Simeonov et al. data hmid barcode file into an indel matrix format.

    This code reads a tab-separated input file containing barcode information, processes the data to
    generate a matrix representation, and writes the results to the specified output file. Additionally,
    it generates two supplementary files: one mapping barcodes to tissues and another mapping mutations.
    """

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    df = pd.read_csv(input_file, sep="\t", index_col=0)
    df.index = ["barcode" + str(idx) for idx in df.index]

    # Use cells per barcode to determine tissues
    tissues = {}
    cell_barcodes = {}
    num_cells = {}
    for idx, row in df.iterrows():
        tissues[idx] = set([cell.strip().split("_")[0] for cell in row["cells"].split(",")])
        cell_barcodes[idx] = [cell.strip() for cell in row["cells"].split(",")]
        num_cells[idx] = len(cell_barcodes[idx])

    df = df.drop(columns=["cells"])

    # Split barcode string into site columns
    hmid_split = df["hmid"].str.split("-", expand=True)
    hmid_split.columns = [f"r{i+1}" for i in range(hmid_split.shape[1])]
    df = df.drop(columns=["hmid"])
    df = pd.concat([df, hmid_split], axis=1)

    # Replace UNKNOWN with -1 and NONE with 0
    df = df.replace("UNKNOWN", -1)
    df = df.replace("NONE", 0)

    # Now replace mut string with integers successively while also assigning multiple site mut str to the first site with other sites as -1 (needed for long deletions)
    mut_mapping = {}
    next_mut_id = 1
    mut_freqs = {}
    for idx, row in df.iterrows():
        prev_val = None
        for col in df.columns:
            val = row[col]
            # Skip unedited or missing data sites
            if val == -1 or val == 0:
                prev_val = None
            # if val == prev_val, set to -1
            elif val == prev_val:
                df.at[idx, col] = -1
            elif val in mut_mapping:
                df.at[idx, col] = mut_mapping[val]
                mut_freqs[mut_mapping[val]] += num_cells[idx]    # Compute frequencies weighted by number of cells with the edit
                prev_val = val
            else:
                mut_mapping[val] = next_mut_id
                df.at[idx, col] = next_mut_id
                mut_freqs[next_mut_id] = num_cells[idx]  # Initialize frequency count by num cells
                next_mut_id += 1
                prev_val = val
    
    # Calculate mut frequencies
    total_edits = sum(mut_freqs.values())
    mut_freqs = {k: v / total_edits for k, v in mut_freqs.items()}
    mut_freqs_only = [freq for mut_val, freq in dict(sorted(mut_freqs.items(), key=lambda item: item[0])).items()]
    
    # Get an expanded mutation matrix with rows for each unique tissue per barcode
    tissues_df = pd.DataFrame({
        'group_name': list(tissues.keys()),
        'tissues': [','.join(sorted(tissue_set)) for tissue_set in tissues.values()]
    })
    expanded_matrix_df, expanded_tissues_df = expand_clones_with_multiple_tissues(df, tissues_df)

    # Set output file paths
    output_tissues_file = output_file.replace("_matrix.csv", "_tissues.txt")
    output_mut_mapping_file = output_file.replace("_matrix.csv", "_mut_mapping.txt")
    output_cell_barcode_file = output_file.replace("_matrix.csv", "_cell_barcodes.txt")
    output_mut_freqs_file = output_file.replace("_matrix.csv", "_mut_freqs.txt")
    output_mut_priors_file = output_file.replace("_matrix.csv", "_mut_priors.txt")
    
    output_matrix_expanded_file = output_file.replace("_matrix.csv", "_matrix_expanded.tsv")
    output_tissues_expanded_file = output_file.replace("_matrix.csv", "_tissues_expanded.txt")
    output_fasta_expanded_file = output_file.replace("_matrix.csv", "_matrix_expanded.fasta")
    output_fasta_expanded_file_beam = output_file.replace("_matrix.csv", "_matrix_expanded_beam.fasta")  # Just for convenience since beam expects -1 to be replaced with the next integer
    
    # Write all to output csv files
    with open(output_mut_mapping_file, "w") as f:
        f.write("successive_char_int,mut_str\n")
        for mut_str, mut_id in mut_mapping.items():
            f.write(f"{mut_id},{mut_str}\n")

    with open(output_tissues_file, "w") as f:
        f.write("group_name,tissues\n")
        for barcode, tissue_set in tissues.items():
            f.write(f"{barcode},{';'.join(sorted(tissue_set))}\n")
    
    with open(output_cell_barcode_file, "w") as f:
        f.write("group_name,cell_barcodes\n")
        for barcode, cell_list in cell_barcodes.items():
            f.write(f"{barcode},{';'.join(sorted(cell_list))}\n")
    
    with open(output_mut_freqs_file, "w") as f:
        for freq in mut_freqs_only:
            f.write(f"{freq:.15f}\n")
    
    with open(output_mut_priors_file, "w") as f:
        f.write("mutation_code,rate\n")
        for key, freq in mut_freqs.items():
            f.write(f"{key},{freq:.15f}\n")

    df.to_csv(output_file, index=True, header=True)
    
    expanded_matrix_df.to_csv(output_matrix_expanded_file, index=True, header=True, sep='\t')
    expanded_tissues_df.to_csv(output_tissues_expanded_file, index=False, header=False)
    
    with open(output_fasta_expanded_file, "w") as f:
        for idx, row in expanded_matrix_df.iterrows():
            f.write(f">{idx}\n")
            seq_str = ','.join([str(val) for val in row.values])
            f.write(f"{seq_str}\n")
    
    largest_edit_id = expanded_matrix_df.values.max()
    expanded_matrix_df_beam = expanded_matrix_df.replace(-1, largest_edit_id + 1)
    with open(output_fasta_expanded_file_beam, "w") as f:
        for idx, row in expanded_matrix_df_beam.iterrows():
            f.write(f">{idx}\n")
            seq_str = ','.join([str(val) for val in row.values])
            f.write(f"{seq_str}\n")

simeonov_main_data_dir = "/grid/siepel/home/staklins/projects/crispr_barcode/data/simeonov_2021_real_data/mendeley_data"
m1_data_dir = simeonov_main_data_dir + "/M1/TreeUtils/clone_hmids"
m2_data_dir = simeonov_main_data_dir + "/M2/TreeUtils/clone_hmids"


outdir = "/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/simeonov_preprocess_2021_pancreatic_cancer_data/reformatted_raw_data/"
os.makedirs(outdir, exist_ok=True)

m1_hmid_files = glob(f"{m1_data_dir}/*.txt")
m2_hmid_files = glob(f"{m2_data_dir}/*.txt")

for hmid_file in m1_hmid_files:
    clone_name = "_".join(os.path.basename(hmid_file).split("_")[:2])
    output_file = f"{outdir}/m1_{clone_name}_matrix.csv"
    convert_simeonov_barcode_file_to_indel_matrix(hmid_file, output_file)

for hmid_file in m2_hmid_files:
    clone_name = "_".join(os.path.basename(hmid_file).split("_")[:2])
    output_file = f"{outdir}/m2_{clone_name}_matrix.csv"
    convert_simeonov_barcode_file_to_indel_matrix(hmid_file, output_file)
    

# Get some summary stats on the clones
all_matrix_files = glob(f"{outdir}/*_matrix.csv")
summary_rows = []
for matrix_file in all_matrix_files:
    df = pd.read_csv(matrix_file, sep=",", index_col=0)
    clone_name = os.path.basename(matrix_file).replace("_matrix.csv", "")
    num_cells = df.shape[0]
    num_sites = df.shape[1]
    num_edits = df.values[df.values > 0].shape[0]
    num_informative = sum(df.apply(count_informative_characters))
    tissue_file = matrix_file.replace("_matrix.csv", "_tissues.txt")
    tissue_df = pd.read_csv(tissue_file, sep=",")
    num_tissues = tissue_df['tissues'].apply(lambda x: len(x.split(';'))).nunique()
    summary_rows.append({'clone_name': clone_name, 
                         'num_cells': num_cells, 
                         'num_sites': num_sites, 
                         'num_edits': num_edits, 
                         'num_tissues': num_tissues,
                         'num_phylogenetically_informative_muts': num_informative})

summary_df = pd.DataFrame(summary_rows)
summary_df.sort_values(by='num_cells', ascending=False, inplace=True)

summary_outdir = os.path.dirname(outdir)
summary_df.to_csv(f"{summary_outdir}/simeonov_clone_summary_stats.csv", index=False, sep='\t')

# Subset to interesting clones that will have phylogenetic information and tissues to model migration graphs for
subset_df = summary_df[
    (summary_df['num_cells'] > 2) & 
    (summary_df['num_cells'] < 300) & 
    (summary_df['num_edits'] > 0) & 
    (summary_df['num_tissues'] > 1) & 
    (summary_df['num_phylogenetically_informative_muts'] > 0)
    ]
subset_df.to_csv(f"{summary_outdir}/simeonov_clone_summary_stats_subset.csv", index=False, sep='\t')

