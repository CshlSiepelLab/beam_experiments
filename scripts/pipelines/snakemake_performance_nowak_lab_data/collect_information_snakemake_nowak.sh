#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam_information/"

outfile="$main_dir/all_mouse_cp_information.csv"

echo "mmus,cp,information" > $outfile

# Assumes only one particle per MMUS/CP combination, otherwise need to first combine chains with NSLogAnalyzer and then run this script on the combined terminal log
files=$(find $main_dir -type f -name "terminal*")

for file in $files; do
    mmus=$(echo $file | rev | cut -d'/' -f4 | rev)
    cp=$(echo $file | rev | cut -d'/' -f3 | rev)
    information=$(grep "Information:" $file | tail -n 1 | cut -d' ' -f2)
    echo "$mmus,$cp,$information" >> $outfile
done

