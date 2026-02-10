
import os
import re
from glob import glob
import pandas as pd
import pickle
import cassiopeia as cas

from beam_sup.matrix_utils import convert_matrix_to_row_successive_matrix, count_informative_characters, expand_clones_with_multiple_tissues, collapse_character_matrix

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
mouse_ids = []
for file in met_char_matrix_files:
    mouse_id = "_".join(os.path.basename(file).split('_')[0:3])
    mouse_ids.append(mouse_id)

# Re-do preprocessing
alleletable_infile = "/grid/siepel/home/staklins/projects/crispr_barcode/data/yang_2022_real_data/yang_2022_cell_KPTracer-Data/KPTracer.alleleTable.FINAL.txt"  # For barcode data
allele_table = pd.read_csv(alleletable_infile, sep="\t", index_col=0)

metadata_file = "/grid/siepel/home/staklins/projects/crispr_barcode/data/yang_2022_real_data/metadata/KPTracer_meta.csv"    # For tissue label and time data
metadata = pd.read_csv(metadata_file, index_col=0)

# Drop nan SubTumor entries from metadata
metadata = metadata.dropna(subset=['SubTumor'])

# Preprocess tissue labels
metadata['SubTumor'] = ["_".join(subtumor.split('_')[2:]) for subtumor in metadata['SubTumor']] # Ignore mouse_genotype prefixes
metadata['SubTumor'] = [subtumor.split('_')[0] for subtumor in metadata['SubTumor'] if "_" in subtumor or subtumor] # Collapse multiple tumor pieces from the same tumor
metadata['SubTumor'] = ['T' if subtumor.startswith('T') else subtumor for subtumor in metadata['SubTumor']] # Collapse multiple primary tumors into one
metadata['SubTumor'] = [re.sub(r'\d+', '', subtumor) for subtumor in metadata['SubTumor']]  # Remove digits from tissue labels (specifically collapses multiple tumors or subtumor pieces in the same tissue to that tissue label), preferring coarser labels for consistency with what I did with Quinn et al data.

char_matrices = {}
mut_dicts = {}
collapsed_char_matrices = {}
collapsing_dicts = {}
expanded_collapsed_char_matrices = {}
expanded_tissues_dfs = {}
tissue_label_dicts = {}
experiment_times = {}
successive_matrices = {}
new_mut_dicts = {}
successive_edit_rates = {}
successive_mut_dicts = {}

for lineage in mouse_ids:
    print(f"\n\nProcessing lineage: {lineage}")
    if "All" in lineage:    # Modified to match all sub-families since there is likely not one metastatic fam deposited for one intBC, per Matt's explanation of "Fam" vs "All"
        all_lineage = "_".join(lineage.split("_")[0:2])[:-1]    # Remove the _T# at the end to keep all primaries
        group = allele_table[allele_table["MetFamily"].str[:len(all_lineage)] == all_lineage]
    else:
        group = allele_table[allele_table["MetFamily"] == lineage]
    indel_priors = cas.pp.compute_empirical_indel_priors(group)
    char_matrix_df, priors, mut_dicts[lineage] = cas.pp.convert_alleletable_to_character_matrix(
            group,
            missing_data_state="-1",
            allele_rep_thresh=0.98, # 0.98 from Yang et al. methods section "Calling clonal populations and creating character matrices"
            mutation_priors=indel_priors,
        )
    # Remove columns with > 80% of entries as missing, again from Yang et al. methods
    columns_to_keep = char_matrix_df.columns[char_matrix_df.isin(['-1']).mean() <= 0.8]
    char_matrix_df = char_matrix_df[columns_to_keep]
    char_matrices[lineage] = char_matrix_df
    
    # Get tissues for each cell in the character matrix
    metadata_subset = metadata.loc[char_matrix_df.index.tolist()]
    assert metadata_subset.shape[0] == char_matrix_df.shape[0], f"Metadata subset size does not match character matrix for family {lineage}, check if cells were lost."
    times = metadata_subset['Aging_Time'].unique()
    assert len(times) == 1, f"Multiple aging times found for family {lineage}, which should not happen."
    experiment_times[lineage] = float(times[0])
    tissue_label_dicts[lineage] = metadata_subset['SubTumor'].to_dict()
    
    # Collapse clones with identical character states but different tissue labels
    collapsed_char_matrices[lineage], collapsing_dicts[lineage] = collapse_character_matrix(char_matrix_df, tissue_label_dicts[lineage])
    
    # Expand clones with multiple tissues into multiple rows, one per tissue
    tissue_label_df = pd.DataFrame.from_dict(collapsing_dicts[lineage]['group_to_tissues'], orient='index', columns=['tissues'])
    tissue_label_df.index.name = 'group_name'
    tissue_label_df = tissue_label_df.reset_index()
    expanded_collapsed_char_matrices[lineage], expanded_tissues_dfs[lineage] = expand_clones_with_multiple_tissues(collapsed_char_matrices[lineage], tissue_label_df)
    
    # Make matrix successive
    successive_matrices[lineage], successive_mut_dicts[lineage], successive_edit_rates[lineage] = convert_matrix_to_row_successive_matrix(expanded_collapsed_char_matrices[lineage], mut_dicts[lineage], indel_priors)

# Output raw character matrices
for fam, char_matrix in char_matrices.items():
    char_matrix.to_csv(f"{outdir}/{fam}_raw_character_matrix.tsv", sep="\t")

# Output original mut dict
for fam, mut_dict in mut_dicts.items():
    with open(f"{outdir}/{fam}_raw_mutation_dict.tsv", "w") as f:
        f.write("site,int,mut_str\n")
        for site, muts in mut_dict.items():
            for char_int, mut_str in muts.items():
                f.write(f"{site},{char_int},{mut_str}\n")

# Output collapsed character matrices
for fam, collapsed_char_matrix in collapsed_char_matrices.items():
    collapsed_char_matrix.to_csv(f"{outdir}/{fam}_collapsed_character_matrix.tsv", sep="\t")

# Output collapsing dict tissues
for fam, collapsing_dict in collapsing_dicts.items():
    with open(f"{outdir}/{fam}_collapsed_tissues.tsv", "w") as f:
        for clone, tissues in collapsing_dict['group_to_tissues'].items():
            f.write(f"{clone}\t{tissues}\n")

# Output collapsing dict originals
for fam, collapsing_dict in collapsing_dicts.items():
    with open(f"{outdir}/{fam}_collapsed_originals.tsv", "w") as f:
        for clone, originals in collapsing_dict['group_to_originals'].items():
            f.write(f"{clone}\t{originals}\n")

# Output expanded collapsed character matrices
for fam, expanded_collapsed_char_matrix in expanded_collapsed_char_matrices.items():
    expanded_collapsed_char_matrix.to_csv(f"{outdir}/{fam}_expanded_collapsed_character_matrix.tsv", sep="\t")

# Output expanded tissues dfs
for fam, expanded_tissues_df in expanded_tissues_dfs.items():
    expanded_tissues_df.to_csv(f"{outdir}/{fam}_expanded_collapsed_tissues.csv", sep=",", index=False, header=False)

# Output successive matrices
for fam, successive_matrix in successive_matrices.items():
    successive_matrix.to_csv(f"{outdir}/{fam}_successive_character_matrix.tsv", sep="\t")

# Output successive matrix in laml format
for fam, successive_matrix in successive_matrices.items():
    successive_matrix_laml = successive_matrix.astype(str).replace("-1", "?")
    successive_matrix_laml.to_csv(f"{outdir}/{fam}_successive_character_matrix_laml.csv")

# Output successive matrix in fasta format for beam
for fam, successive_matrix in successive_matrices.items():
    max_mut_val = successive_matrix.astype(int).values.flatten().max()
    new_missing_state = max_mut_val + 1
    successive_matrix_beam = successive_matrix.astype(str).replace("-1", f"{str(new_missing_state)}")
    with open(f"{outdir}/{fam}.fasta", "w") as f:
        for idx, row in successive_matrix_beam.iterrows():
            f.write(f">{idx}\n")
            row_str = ",".join([str(state) for state in row.tolist()])
            f.write(f"{row_str}\n")

# Output successive edit rates
for fam, edit_rates in successive_edit_rates.items():
    sorted_edit_rates = dict(sorted(edit_rates.items(), key=lambda item: item[0]))
    with open(f"{outdir}/{fam}_successive_mut_priors.txt", "w") as f:
        f.write("mutation_code,rate\n")
        for mut, rate in sorted_edit_rates.items():
            f.write(f"{mut},{rate}\n")
    just_sorted_rates = [rate for mut, rate in sorted_edit_rates.items()]
    with open(f"{outdir}/{fam}_successive_mut_freqs.txt", "w") as f:
        for rate in just_sorted_rates:
            f.write(f"{rate}\n")

# Output successive mut dict
for fam, mut_dict in successive_mut_dicts.items():
    with open(f"{outdir}/{fam}_successive_mutation_dict.tsv", "w") as f:
        f.write("successive_char_int,mut_str\n")
        for mut_str, successive_char_int in mut_dict.items():
            f.write(f"{successive_char_int},{mut_str}\n")

# Output experiment times
with open(f"{outdir}/met_fam_experiment_times.txt", "w") as f:
    f.write("met_fam_id\texperiment_time\n")
    for fam, time in experiment_times.items():
        f.write(f"{fam}\t{time}\n")

# Get number of cells per family in the raw data (not collapsed)
num_cells_per_family = {fam: df.shape[0] for fam, df in char_matrices.items()}
with open(f"{outdir}/met_fam_num_cells_raw.txt", "w") as f:
    f.write("met_fam_id\tnum_cells\n")
    for fam, num_cells in num_cells_per_family.items():
        f.write(f"{fam}\t{num_cells}\n")

num_cells_per_family_collapsed = {fam: df.shape[0] for fam, df in collapsed_char_matrices.items()}
with open(f"{outdir}/met_fam_num_cells_collapsed.txt", "w") as f:
    f.write("met_fam_id\tnum_cells\n")
    for fam, num_cells in num_cells_per_family_collapsed.items():
        f.write(f"{fam}\t{num_cells}\n")


unique_tissues_per_family = {fam: set(tissues.values()) for fam, tissues in tissue_label_dicts.items()}
with open(f"{outdir}/met_fam_unique_tissue_labels.txt", "w") as f:
    f.write("met_fam_id\ttissues\n")
    for fam, tissues in unique_tissues_per_family.items():
        f.write(f"{fam}\t{','.join(tissues)}\n")
            


# Get some summary stats on the clones
summary_rows = []
for fam, df in collapsed_char_matrices.items():
    df = df.astype(int)
    num_cells = df.shape[0]
    num_sites = df.shape[1]
    num_edits = df.values[df.values > 0].shape[0]
    num_informative = sum(df.apply(count_informative_characters))
    num_tissues = expanded_tissues_dfs[fam]['tissues'].nunique()
    summary_rows.append({'clone_name': fam, 
                         'num_cells': num_cells, 
                         'num_sites': num_sites, 
                         'num_edits': num_edits, 
                         'num_tissues': num_tissues,
                         'num_phylogenetically_informative_muts': num_informative})

summary_df = pd.DataFrame(summary_rows)
summary_df.sort_values(by='num_cells', ascending=False, inplace=True)
summary_df.to_csv(f"{outdir}/yang_clone_summary_stats.tsv", index=False, sep='\t')

# Subset to interesting clones that will have phylogenetic information and tissues to model migration graphs for
subset_df = summary_df[
    (summary_df['num_cells'] > 2) & 
    (summary_df['num_cells'] < 300) & 
    (summary_df['num_edits'] > 0) & 
    (summary_df['num_tissues'] > 1) & 
    (summary_df['num_phylogenetically_informative_muts'] > 0)
    ]
subset_df.to_csv(f"{outdir}/yang_clone_summary_stats_subset.tsv", index=False, sep='\t')


