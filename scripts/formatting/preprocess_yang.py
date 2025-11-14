
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
# with open(f"{outdir}/met_fam_num_cells_raw.txt", "w") as f:
#     f.write("met_fam_id\tnum_cells\n")
#     for fam, num_cells in num_cells_per_family.items():
#         f.write(f"{fam}\t{num_cells}\n")

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
# with open(f"{outdir}/met_fam_unique_tissue_labels.txt", "w") as f:
#     f.write("met_fam_id\ttissues\n")
#     for fam, tissues in unique_tissues_per_family.items():
#         f.write(f"{fam}\t{','.join(tissues)}\n")

# Collapse the character matrices to unique rows and create collapsing dicts for those clones
collapsed_char_matrices = {}
collapsing_dicts = {}
for fam, df in char_matrices.items():
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
        original_tissues = ",".join(list(set([tissue_label_dicts[fam][cell] for cell in original_row_names])))
        group_to_tissues[group_name] = original_tissues
        
    # Replace index names in unique_rows with the appropriate group name
    unique_rows.index = group_names
    unique_rows = unique_rows.replace('-', -1)  # Ensure missing data is -1 as integer
    collapsed_char_matrices[fam] = unique_rows
    collapsing_dicts[fam] = {'group_to_originals': group_to_originals, 'group_to_tissues': group_to_tissues}

num_cells_per_family_collapsed = {fam: df.shape[0] for fam, df in collapsed_char_matrices.items()}
# with open(f"{outdir}/met_fam_num_cells_collapsed.txt", "w") as f:
#     f.write("met_fam_id\tnum_cells\n")
#     for fam, num_cells in num_cells_per_family_collapsed.items():
#         f.write(f"{fam}\t{num_cells}\n")

# # Read in matched priors files from deposited data by cassiopeia analysis
# priors_files = [file.replace("_character_matrix.txt", "_priors.pkl") for file in met_char_matrix_files]
# priors_dicts = {}
# for file in priors_files:
#     mouse_id = "_".join(os.path.basename(file).split('_')[0:3])
#     with open(file, "rb") as f:
#         priors_dicts[mouse_id] = pickle.load(f)

# assert priors_dicts.keys() == char_matrices.keys(), "Priors and character matrices do not match in metastatic families"

# # Make all site names start at r1, r2, ...
# # This is just my preference for what I've been using elsewhere with quinn et al. and simeonov et al. data
# for fam, df in collapsed_char_matrices.items():
#     df.columns = [f"r{int(col[1:]) + 1}" for col in df.columns]
#     priors_dicts[fam] = {f"r{key + 1}": vals for key, vals in priors_dicts[fam].items()}
#     assert df.columns.tolist() == list(priors_dicts[fam].keys()), f"Column names do not match priors keys for family {fam}"

# def get_successive_matrix(character_matrix, priors):
#     """ New function here due to slightly different priors format in the deposited data vs my previous function in beam_sup package """
#     successive_char_matrix = character_matrix.copy()
#     successive_mut_dict = {}
#     i = 1
#     for clone, row in character_matrix.iterrows():
#         for site, mut in row.items():
#             mut = int(mut)
#             # Skip undedited and missing sites
#             if mut == 0 or mut == -1:
#                 continue
#             mut = str(mut)
#             # Replace the mutation with the successive mutation
#             if mut not in successive_mut_dict:
#                 successive_mut_dict[mut] = i
#                 new_mut_value = i
#                 i += 1
#             else:
#                 new_mut_value = successive_mut_dict[mut]
#             successive_char_matrix.loc[clone, site] = new_mut_value
            
#     # For maintaining laml style priors
#     successive_priors = {}
#     for site, muts in priors.items():
#         successive_priors[site] = {}
#         for mut, mut_prior_freq in muts.items():
#             mut_int = successive_mut_dict[str(mut)]
#             successive_priors[site][mut_int] = mut_prior_freq
            
#     return successive_char_matrix, successive_mut_dict, successive_priors

# # Convert collapsed char matrices to successive format and update priors accordingly
# successive_char_matrices = {}
# successive_priors = {}
# successive_mut_mappings = {}
# for fam, df in collapsed_char_matrices.items():
#     successive_char_matrices[fam], successive_mut_mappings[fam], successive_priors[fam] = get_successive_matrix(df, priors_dicts[fam])
    
# # Normalize successive priors site-wise since cassiopeia does not do this automatically but laml expects probabilities per site implying a sum to 1
# successive_priors_normalized = {}
# for fam, priors in successive_priors.items():
#     successive_priors_normalized[fam] = {}
#     for site, muts in priors.items():
#         total_rate = sum(muts.values())
#         successive_priors_normalized[fam][site] = {mut: rate / total_rate for mut, rate in muts.items()}

# Get edit rates for beam across the full matrix all sites for the sequential frequencies
# # Option 1: get frequencies from the observed edits in the full matrix weighted by number of cells in each clone row
# full_matrix_edit_freqs = {}
# for fam, df in successive_char_matrices.items():
#     full_matrix_edit_freqs[fam] = {}
#     for clone, row in df.iterrows():
#         clone_num_cells = len(collapsing_dicts[fam]['group_to_originals'][clone].split(','))
#         for site, mut in row.items():
#             if mut == -1 or mut == 0:
#                 continue
#             if mut not in full_matrix_edit_freqs[fam]:
#                 full_matrix_edit_freqs[fam][mut] = clone_num_cells
#             else:
#                 full_matrix_edit_freqs[fam][mut] += clone_num_cells    # Keep track of the number of barcodes with this edit
#     # Convert counts to frequencies
#     total_num_muts = sum(full_matrix_edit_freqs[fam].values())
#     full_matrix_edit_freqs[fam] = {mut: count / total_num_muts for mut, count in full_matrix_edit_freqs[fam].items()}

# # Option 2: repurpose the edit rates from the deposited data priors
# full_matrix_repurposed_edit_rates = {}
# for fam, priors in successive_priors.items():
#     full_matrix_repurposed_edit_rates[fam] = {}
#     for site, muts in priors.items():
#         for mut, rate in muts.items():
#             if mut not in full_matrix_repurposed_edit_rates[fam]:
#                 full_matrix_repurposed_edit_rates[fam][mut] = rate
#             else:
#                 full_matrix_repurposed_edit_rates[fam][mut] += rate  # Sum rates across sites
#     # Normalize to sum to 1
#     total_rate = sum(full_matrix_repurposed_edit_rates[fam].values())
#     full_matrix_repurposed_edit_rates[fam] = {mut: rate / total_rate for mut, rate in full_matrix_repurposed_edit_rates[fam].items()}
            

# Now I will try to re-do the preprocessing since the deposited data above does not have the mutation info across sites
alleletable_infile = "/grid/siepel/home/staklins/projects/crispr_barcode/data/yang_2022_real_data/yang_2022_cell_KPTracer-Data/KPTracer.alleleTable.FINAL.txt"

allele_table = pd.read_csv(alleletable_infile, sep="\t")

groups_from_deposited_data = list(char_matrices.keys())

for lineage in groups_from_deposited_data:
        if "All" in lineage:
            continue    # Just skip for now
        group = allele_table[allele_table["MetFamily"] == lineage]

        indel_priors = cas.pp.compute_empirical_indel_priors(group)

        char_matrix_df, priors, mut_dict = cas.pp.convert_alleletable_to_character_matrix(
                group,
                missing_data_state="-1",
                allele_rep_thresh=0.98, # From Yang et al. methods section "Calling clonal populations and creating character matrices"
                mutation_priors=indel_priors,
            )
        assert char_matrix_df.shape == char_matrices[lineage].shape, "Character matrix shape does not match deposited data for lineage {lineage}"

        # Remove columns with > 80% of entries as missing, again from Yang et al. methods
        columns_to_keep = char_matrix_df.columns[char_matrix_df.isin(['-1']).mean() <= 0.8]
        char_matrix_df = char_matrix_df[columns_to_keep]

        successive_matrix, new_mut_dict, successive_edit_rates = convert_matrix_to_row_successive_matrix(char_matrix_df, mut_dict, indel_priors)

### Old code below whre I did my own preprocessing, but the resulting matrices were different so I reverted to repurposing the existing ones from the deposited data above
# os.makedirs(outdir, exist_ok=True)

# allele_df = pd.read_csv(allele_filepath, sep="\t", index_col=0)
# # Find metastatic mice
# # Make new column with tissue labels only
# allele_df["tissue"] = [
# "".join(filter(str.isalpha, name[2]))
# for name in allele_df["Tumor"].str.split("_")
# ]

# # Get met lineage names and tissues
# tumor_names = allele_df.groupby("MetFamily")["tissue"].unique()

# # Keep only lineages with both primary and met tissues
# tumor_names = tumor_names[
# tumor_names.apply(
#     lambda x: any(name.startswith("T") for name in x)
#     and any(not name.startswith("T") for name in x)
# )
# ]

# # keep only subset MetFamily mice
# allele_df = allele_df[allele_df["MetFamily"].isin(tumor_names.index)]

# # Make character matrix
# # Write mouse specific character matrix to file
# for MetFamily, df in allele_df.groupby("MetFamily"):

# # get indel priors as per Cassiopeia docs
# indel_priors = cas.pp.compute_empirical_indel_priors(df)

# # Built in cassiopeia function to convert allele table to character matrix
# allele_thresh = 0.98    # From Yang et al. methods
# char_matrix_df, prob_dict, mut_dict = cas.pp.convert_alleletable_to_character_matrix(
#     df,
#     missing_data_state="-1",
#     allele_rep_thresh=allele_thresh,
#     mutation_priors=indel_priors
# )

# successive_matrix, new_mut_dict, successive_edit_rates = (
#     convert_matrix_to_row_successive_matrix(char_matrix_df, mut_dict, indel_priors)
# )

# # DID NOT UPDATE PAST THIS POINT YET FOR NEW SUCCESSIVE MATRIX FUNCTION

# # Rename columns of successive char matrix to be successive themselves
# successive_char_matrix.columns = [
#     f"r{i}" for i in range(1, len(successive_char_matrix.columns) + 1)
# ]

# # Output successive char matrix for all cells
# mouse_outfile = f"{outdir}/{MetFamily}_successive_char_matrix.txt"
# successive_char_matrix.index.name = "cellBC"
# # successive_char_matrix.to_csv(mouse_outfile, sep="\t", index=True)

# # Write mutation dictionary to file
# mut_dict_outfile = mouse_outfile.replace(".txt", f"_mut_dict.txt")
# # with open(mut_dict_outfile, "w") as f:
# #     f.write(f"mut_id\tmut_str\n")
# #     for str, id in successive_mut_dict.items():
# #         f.write(f"{id}\t{str}\n")

# # Collapse the cells to only unique rows and output collapsing dict of cellBCs and tissue labels
# all_columns = successive_char_matrix.columns.tolist()
# sorted_char_matrix = successive_char_matrix.sort_values(by=all_columns)
# unique_rows = sorted_char_matrix.drop_duplicates(keep="first")
# group_names = [f"clone{i+1}" for i in range(len(unique_rows))]
# group_to_originals = {}
# group_to_tissues = {}
# for group_name, (_, unique_row) in zip(group_names, unique_rows.iterrows()):
#     # Find all rows in sorted_char_matrix that match the unique_row
#     original_row_names = sorted_char_matrix[
#         sorted_char_matrix.eq(unique_row).all(axis=1)
#     ].index.tolist()
#     group_to_originals[group_name] = original_row_names
#     original_tissues = set(
#         df[df["cellBC"].isin(original_row_names)]["tissue"].values.tolist()
#     )
#     group_to_tissues[group_name] = original_tissues

# # Replace index names in unique_rows with the appropriate group name
# unique_rows.index = group_names

# print(
#     MetFamily,
#     f"cells: {len(successive_char_matrix)}",
#     f"clones: {len(unique_rows)}",
#     f"sites: {len(successive_char_matrix.columns)}",
# )

# # Write unique rows to file
# unique_rows_outfile = mouse_outfile.replace(".txt", f"_collapsed.txt")
# # unique_rows.to_csv(unique_rows_outfile, sep="\t", index=True)

# # Write collapsing dict of cellBCs and tissue labels to file
# collapsing_dict_outfile = mouse_outfile.replace(".txt", f"_collapsing_dict.txt")
# # with open(collapsing_dict_outfile, "w") as f:
# #     f.write(f"group_name\tcellBCs\ttissues\n")
# #     for group_name in group_names:
# #         cellBCs = ','.join(list(group_to_originals[group_name]))
# #         tissues = ','.join(list(group_to_tissues[group_name]))
# #         f.write(f"{group_name}\t{cellBCs}\t{tissues}\n")
