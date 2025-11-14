
import os
from glob import glob
import pandas as pd
import pickle
import cassiopeia as cas

from beam_sup.matrix_utils import convert_matrix_to_row_successive_matrix, count_informative_characters, expand_clones_with_multiple_tissues

# Where to output processed files
outdir = "/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/yang_preprocess_2022_lung_cancer_data"
os.makedirs(outdir, exist_ok=True)

# Find the metastatic met family character matrices from the deposited data
deposited_char_matrix_path="/grid/siepel/home/staklins/projects/crispr_barcode/data/yang_2022_real_data/yang_2022_cell_KPTracer-Data/trees"
char_matrix_files = glob(f"{deposited_char_matrix_path}/*_character_matrix.txt")
met_char_matrix_files = [f for f in char_matrix_files if "all" in os.path.basename(f).lower() or "fam" in os.path.basename(f).lower()]

# Filter to remove duplicate where both All and Fam versions exist, preferring the smaller Fam versions
fam_files = [f for f in met_char_matrix_files if "fam" in os.path.basename(f).lower()]
removed_duplicates = []
for f in fam_files:
    file_prefix = "_".join(f.split("_")[:-4])  # Remove the last four parts (_T#_Fam_character_matrix.txt)
    all_file = file_prefix + "_All_character_matrix.txt"
    if all_file in met_char_matrix_files:
        met_char_matrix_files.remove(all_file)
        removed_duplicates.append(all_file)
    
# Read in the character matrices for metastatic families
char_matrices = {}
for file in met_char_matrix_files:
    mouse_id = "_".join(os.path.basename(file).split('_')[0:3])
    char_matrices[mouse_id] = pd.read_csv(file, sep='\t', index_col=0)

# Get number of cells per family in the raw data (not collapsed)
num_cells_per_family = {fam: df.shape[0] for fam, df in char_matrices.items()}
with open(f"{outdir}/met_fam_num_cells_raw.txt", "w") as f:
    f.write("met_fam_id\tnum_cells\n")
    for fam, num_cells in num_cells_per_family.items():
        f.write(f"{fam}\t{num_cells}\n")

# Get the tissue labels for each cell from the metadata file for metastatic families
metadata_file = "/grid/siepel/home/staklins/projects/crispr_barcode/data/yang_2022_real_data/metadata/KPTracer_meta.csv"
metadata = pd.read_csv(metadata_file, index_col=0)

# Get the tissue labels from the SubTumor sample names
all_cell_names = []
for fam, df in char_matrices.items():
    all_cell_names.extend(df.index.tolist())
    
metadata = metadata.loc[all_cell_names]
assert metadata.shape[0] == len(all_cell_names), "Metadata does not match character matrix cell names, check if cells were lost."

# Preprocess tissue labels
metadata['SubTumor'] = ["_".join(subtumor.split('_')[2:]) for subtumor in metadata['SubTumor']] # Ignore mouse_genotype prefixes
metadata['SubTumor'] = [subtumor.split('_')[0] for subtumor in metadata['SubTumor'] if "_" in subtumor or subtumor] # Collapse multiple tumor pieces from the same tumor
metadata['SubTumor'] = ['T' if subtumor.startswith('T') else subtumor for subtumor in metadata['SubTumor']] # Collapse multiple primary tumors into one

# Get tissue label dicts for each metastatic family
tissue_label_dicts = {}
experiment_times = {}
for fam, df in char_matrices.items():
    metadata_subset = metadata.loc[df.index.tolist()]
    assert metadata_subset.shape[0] == df.shape[0], f"Metadata subset size does not match character matrix for family {fam}"
    times = metadata_subset['Aging_Time'].unique()
    assert len(times) == 1, f"Multiple aging times found for family {fam}"
    experiment_times[fam] = float(times[0])
    tissue_label_dicts[fam] = metadata_subset['SubTumor'].to_dict()

unique_tissues_per_family = {fam: set(tissues.values()) for fam, tissues in tissue_label_dicts.items()}
with open(f"{outdir}/met_fam_unique_tissue_labels.txt", "w") as f:
    f.write("met_fam_id\ttissues\n")
    for fam, tissues in unique_tissues_per_family.items():
        f.write(f"{fam}\t{','.join(tissues)}\n")

def collapse_character_matrices(input_char_matrices, input_tissue_label_dicts):
    """ Collapse character matrices to unique rows and create collapsing dicts for those clones """
    collapsed_char_matrices = {}
    collapsing_dicts = {}
    for fam, df in input_char_matrices.items():
        all_columns = df.columns.tolist()
        sorted_char_matrix = df.sort_values(by=all_columns)
        unique_rows = sorted_char_matrix.drop_duplicates(keep="first")
        group_names = [f"clone{i+1}" for i in range(len(unique_rows))]
        group_to_originals = {}
        group_to_tissues = {}
        for group_name, (_, unique_row) in zip(group_names, unique_rows.iterrows()):
            # Find all rows in sorted_char_matrix that match the unique_row
            original_row_names = sorted_char_matrix[sorted_char_matrix.eq(unique_row).all(axis=1)].index.tolist()
            group_to_originals[group_name] = ",".join(original_row_names)
            original_tissues = ",".join(list(set([input_tissue_label_dicts[fam][cell] for cell in original_row_names])))
            group_to_tissues[group_name] = original_tissues
            
        # Replace index names in unique_rows with the appropriate group name
        unique_rows.index = group_names
        unique_rows = unique_rows.replace('-', -1)  # Ensure missing data is -1 as integer
        collapsed_char_matrices[fam] = unique_rows
        collapsing_dicts[fam] = {'group_to_originals': group_to_originals, 'group_to_tissues': group_to_tissues}
    return collapsed_char_matrices, collapsing_dicts

# Collapse character matrices to unique rows and get collapsing dicts
collapsed_char_matrices, collapsing_dicts = collapse_character_matrices(char_matrices, tissue_label_dicts)

num_cells_per_family_collapsed = {fam: df.shape[0] for fam, df in collapsed_char_matrices.items()}
with open(f"{outdir}/met_fam_num_cells_collapsed.txt", "w") as f:
    f.write("met_fam_id\tnum_cells\n")
    for fam, num_cells in num_cells_per_family_collapsed.items():
        f.write(f"{fam}\t{num_cells}\n")
            

# Now I will try to re-do the preprocessing since the deposited data above does not have the mutation info across sites that I want
alleletable_infile = "/grid/siepel/home/staklins/projects/crispr_barcode/data/yang_2022_real_data/yang_2022_cell_KPTracer-Data/KPTracer.alleleTable.FINAL.txt"

allele_table = pd.read_csv(alleletable_infile, sep="\t", index_col=0)

groups_from_deposited_data = list(char_matrices.keys())

for lineage in groups_from_deposited_data:
    print(f"\n\nProcessing lineage: {lineage}")
    if "All" in lineage:    # Modified to match all sub-families since there is likely not one metastatic fam deposited for one intBC, per Matt's explanation of "Fam" vs "All"
        all_lineage = "_".join(lineage.split("_")[0:2])[:-1]    # Remove the _T# at the end to keep all primaries
        group = allele_table[allele_table["MetFamily"].str[:len(all_lineage)] == all_lineage]
    else:
        group = allele_table[allele_table["MetFamily"] == lineage]
    threshold = 0.98 # 0.98 from Yang et al. methods section "Calling clonal populations and creating character matrices"
    # Manually change some thresholds to match the deposited data above
    if lineage == "3460_Lkb1_T1":
        threshold = 1.0
    if lineage == "3508_Apc_T2":
        threshold = 0.9835
    indel_priors = cas.pp.compute_empirical_indel_priors(group)
    char_matrix_df, priors, mut_dict = cas.pp.convert_alleletable_to_character_matrix(
            group,
            missing_data_state="-1",
            allele_rep_thresh=threshold, 
            mutation_priors=indel_priors,
        )
    # Remove columns with > 80% of entries as missing, again from Yang et al. methods
    columns_to_keep = char_matrix_df.columns[char_matrix_df.isin(['-1']).mean() <= 0.8]
    char_matrix_df = char_matrix_df[columns_to_keep]
    
    if char_matrix_df.shape != char_matrices[lineage].shape:
        print(f"Character matrix shape does not match deposited data for lineage {lineage}")
        print(f"Deposited data shape: {char_matrices[lineage].shape}, Recomputed shape: {char_matrix_df.shape}")
    else:
        print(f"Character matrix shape matches deposited data for lineage {lineage}")
    


successive_matrix, new_mut_dict, successive_edit_rates = convert_matrix_to_row_successive_matrix(char_matrix_df, mut_dict, indel_priors)
