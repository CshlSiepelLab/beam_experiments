#!/bin/bash

dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_pairwise_test_100_tips_9_12_24/precision_recall"

outfile=${dir}/downsample_all_stats_formatted.csv
echo "downsample_threshold,posterior_threshold,precision,recall,sim,thresh_counts" > $outfile

files=$(find $dir -type f -name all_threshold_stats.csv)

for file in $files; do
    threshold=$(dirname $file | rev | cut -d'/' -f1 | rev)
    grep "^0.4," $file | sed "s/^/${threshold},/" >> $outfile
done