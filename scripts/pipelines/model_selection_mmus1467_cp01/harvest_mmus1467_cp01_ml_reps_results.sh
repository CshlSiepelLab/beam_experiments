#!/bin/bash

files=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/beam_ns_mmus1467_cp01_no_reseeding_one_rate_reseeding_12_3_24 -type f -name "*terminal.log")

outfile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/beam_ns_mmus1467_cp01_no_reseeding_one_rate_reseeding_12_3_24/all_reps_results_ml_mmus1467_cp01.csv"

echo -e "model,rep,marginal_likelihood" >> $outfile

for file in $files; do
    echo $file
    name=$(echo $file | rev | cut -d'/' -f2 | rev)
    rep=$(echo $name | rev | cut -d'_' -f1 | rev | sed 's/rep//')
    model=$(echo $name | rev | cut -d'_' -f2- | rev)
    ml=$(grep "Marginal likelihood:" $file | awk -F' ' '{print $3}' | awk -F'(' '{print $1}' | tail -n 1)
    
    echo -e "$model,$rep,$ml" >> $outfile
done