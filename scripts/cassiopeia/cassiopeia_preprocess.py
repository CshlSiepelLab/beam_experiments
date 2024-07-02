#!/usr/bin/env python3

import sys, os
import pandas as pd
import cassiopeia as cas

# USER INPUTS
 # The raw paired FASTQs from a single sample
# input_files = [
#     "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/raw_fastqs/SRR17885786_1.fastq.gz", 
#     "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/raw_fastqs/SRR17885786_2.fastq.gz"
# ]
input_files = sys.argv[1].split(',')
output_dir = sys.argv[2]

# The sample name, used for naming output files
name = os.path.basename(input_files[0]).split('_')[0]

# Path to the target site reference sequence in FASTA format
reference_filepath = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/inputs/PCT48.ref.fasta"

# describe the locations of intBC anad cutsites in the ref sequence
barcode_interval = (20, 34)
cutsite_locations = [112, 166, 220]

# Number of threads to use, whenever parallelization is possible
n_threads = 10

# Whether to allow a single intBC to have multiple allele states. For chemistries for which barcode == cell, this should be `False`.
allow_allele_conflicts = False

# Verbosity of logging
verbose = True

# Specify the version of 10x chemistry used
chem = '10xv2'

# Take from the cellranger default list for the 10x chemistry used
cellbc_whitelist = '/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/inputs/737K-august-2016.txt'

# intbc_whitelist = ''


# MAIN PIPELINE
cas.pp.setup(output_dir, verbose=verbose)

bam_fp = cas.pp.convert_fastqs_to_unmapped_bam(
    input_files,
    chemistry=chem,
    output_directory=output_dir,
    name=name,
    n_threads=n_threads
)
print("Done converting fastqs to unmapped bam")

bam_fp = cas.pp.filter_bam(
    bam_fp,
    output_directory=output_dir,
    quality_threshold=10,
    n_threads=n_threads,
)
print("Done filtering bam")

bam_fp = cas.pp.error_correct_cellbcs_to_whitelist(
    bam_fp,
    whitelist=cellbc_whitelist,
    output_directory=output_dir,
    n_threads=n_threads,
)
print("Done error correcting cell barcodes")

umi_table = cas.pp.collapse_umis(
    bam_fp,
    output_directory=output_dir,
    max_hq_mismatches=3,
    max_indels=2,
    method='likelihood',
    n_threads=n_threads,
)
print("Done collapsing UMIs")

umi_table = cas.pp.resolve_umi_sequence(
    umi_table,
    output_directory=output_dir,
    min_umi_per_cell=10,
    min_avg_reads_per_umi=2.0,
    plot=True,
)
print("Done resolving UMI sequences")

umi_table = cas.pp.align_sequences(
    umi_table,
    ref_filepath=reference_filepath,
    gap_open_penalty=20,
    gap_extend_penalty=1,
    n_threads=n_threads,
)
print("Done aligning sequences")

umi_table = cas.pp.call_alleles(
    umi_table,
    ref_filepath=reference_filepath,
    barcode_interval=barcode_interval,
    cutsite_locations=cutsite_locations,
    cutsite_width=12,
    context=True,
    context_size=5,
)
print("Done calling alleles")

# umi_table = cas.pp.error_correct_intbcs_to_whitelist(
#     umi_table,
#     whitelist=intbc_whitelist,
#     intbc_dist_thresh=1
# )
# print("Done error correcting int barcodes")

umi_table = cas.pp.error_correct_umis(
    umi_table,
    max_umi_distance=2,
    allow_allele_conflicts=allow_allele_conflicts,
    n_threads=n_threads,
)
print("Done error correcting UMIs")

umi_table = cas.pp.filter_molecule_table(
    umi_table,
    output_directory=output_dir,
    min_umi_per_cell=10,
    min_avg_reads_per_umi=2.0,
    min_reads_per_umi=-1,
    intbc_prop_thresh=0.5,
    intbc_umi_thresh=10,
    intbc_dist_thresh=1,
    doublet_threshold=0.35,
    allow_allele_conflicts=allow_allele_conflicts,
    plot=True,
)
print("Done filtering molecule table")

allele_table = cas.pp.call_lineage_groups(
    umi_table,
    output_directory=output_dir,
    min_umi_per_cell=10,
    min_avg_reads_per_umi=2.0,
    min_cluster_prop=0.005,
    min_intbc_thresh=0.05,
    inter_doublet_threshold=0.35,
    kinship_thresh=0.25,
    plot=True,
)
print("Done calling lineage groups")
