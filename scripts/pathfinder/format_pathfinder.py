#!/usr/bin/env python

import sys
import pandas as pd

# user input
tissue_labels = sys.argv[1]
indel_matrix = sys.argv[2]
outdir = sys.argv[3]

# # testing
# tissue_labels = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/raw_data/2945/cell_tree_seed1082116693.labeling'
# indel_matrix = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/raw_data/2945/2945_indel_character_matrix.tsv'
# outdir = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/test/pathfinder/2945'

# read in the data to pd df
df = pd.read_csv(tissue_labels, sep=' ', header=None, names=['Sample', 'Label'])

unique_labels = df['Label'].unique()

# Create a new DataFrame with columns as sample names and rows as site labels
one_hot_df = pd.DataFrame(0, index=unique_labels, columns=df['Sample'].unique())

for _, row in df.iterrows():
    one_hot_df.at[row['Label'], row['Sample']] = 1

one_hot_df.index.name="Tumor"

# Output DataFrame to TSV file
output_file = f"{outdir}/clone_presence.txt"
one_hot_df.to_csv(output_file, sep='\t', index=True)


# format fasta file
indel_df = pd.read_csv(indel_matrix, sep='\t', index_col=0)

# mega does not model missing data, so we need to replace -1 with 0
indel_df.replace(-1, 0, inplace=True)

# get all clone names
site_names = indel_df.columns.values.tolist()
clone_names = indel_df.index.values.tolist()

# get the number of mutations per site to format T/A string
site_sizes = {site: len(indel_df[site].unique()) for site in indel_df.columns}
site_ordered_values = {site: indel_df[site].unique().tolist() for site in indel_df.columns}

total_sequence_size = sum(list(site_sizes.values()))

clone_sequences = {clone: "" for clone in clone_names}
clone_sequences['Normal'] = "A" * total_sequence_size

for clone in clone_names:
    for i, site in enumerate(site_names):
        new_seq = "A" * (site_sizes[site])
        mut = indel_df.at[clone, site]
        if mut == 0:
            pass
        else: 
            nt = site_ordered_values[site].index(mut)
            new_seq = new_seq[:nt] + "T" + new_seq[nt+1:]
        clone_sequences[clone] += new_seq

# Output to fasta file
output_file = f"{outdir}/clone_aln.fas"
with open(output_file, 'w') as f:
    for clone, seq in clone_sequences.items():
        f.write(f">{clone}\n")
        f.write(f"{seq}\n")