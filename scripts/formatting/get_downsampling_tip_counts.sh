#!/bin/bash

# desired output file path
outfile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_10_31_24/downsampled_data/tip_counts.csv"

echo -e "sim_name,downsampling_threshold,tip_count" > $outfile

# path to downsampling dirs
downsampling_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_10_31_24/downsampled_data"

# find files is the downsampling dir
files=$(find $downsampling_dir -name "*_tissue_labels.tsv")

for file in $files; do 
    sim_name=$(echo $file | rev | cut -d'/' -f2 | rev)
    threshold=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'_' -f2)
    num_tips=$(cat $file | wc -l)
    echo -e "$sim_name,$threshold,$num_tips" >> $outfile 
done


