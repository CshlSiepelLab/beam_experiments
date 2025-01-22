#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_1_20_25_uniform_50cells_50sites_data_7_24_24" 
gtr_dir=$main_dir/beam_gtr_ns

outfile=$main_dir/marginal_likelihoods.csv

# Write header
echo "name,gtr_ml,gtr_sd,random_ml,random_sd,bf(gtr-random),diff_threshold" > $outfile

files=$(find $gtr_dir -type f -name "terminal*")

for file in $files; do

    gtr_file=$file
    random_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_random_ns/g')

    if [[ ! -f $random_file || ! -f $gtr_file ]]; then
        echo "Missing file: $random_file or $gtr_file"
        continue
    fi

    name=$(echo $file | rev | cut -d'/' -f2 | rev)

    gtr_ml=$(grep "Marginal likelihood" $gtr_file | tail -n 1 | cut -d' ' -f3 | cut -d'(' -f1)
    gtr_sd=$(grep "Marginal likelihood" $gtr_file | tail -n 1 | cut -d' ' -f3 | cut -d'(' -f2 | sed 's/)//g')

    random_ml=$(grep "Marginal likelihood" $random_file | tail -n 1 | cut -d' ' -f3 | cut -d'(' -f1)
    random_sd=$(grep "Marginal likelihood" $random_file | tail -n 1 | cut -d' ' -f3 | cut -d'(' -f2 | sed 's/)//g')

    # Bayes factor is reported with Hnull as the no reseeding and Halt as one rate reseeding, so a positive Bayes factor value supports reseeding and negative supports no reseeding
    beam_bf=$(echo "scale=10; $gtr_ml - $random_ml" | bc -l)

    # Required Bayes factor difference threshold based on the estimated standard deviations of the marginal likelihoods from nested sampling
    diff_threshold=$(echo "scale=10; 2 * sqrt(($gtr_sd^2) + ($random_sd^2))" | bc -l)

    # Write to outfile
    echo -e "$name,$gtr_ml,$gtr_sd,$random_ml,$random_sd,$beam_bf,$diff_threshold" >> $outfile

done