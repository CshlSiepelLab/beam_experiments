#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/beam_ns_mmus1544_cp04_no_reseeding_one_rate_reseeding_12_18_24"

# No reseeding
dirs=$(find $main_dir -mindepth 1 -maxdepth 1 -type d -name "no_reseeding*")

files=""

for dir in $dirs; do
    name=$(echo $dir | rev | cut -d'/' -f1 | rev)
    files+=" $dir/${name}.log"
done

applauncher NSLogAnalyser -N 1 -noposterior $files -out "$main_dir/no_reseeding_combined.log" > "$main_dir/no_reseeding_combined_terminal.log" 2>&1

# One rate reseeding
dirs=$(find $main_dir -mindepth 1 -maxdepth 1 -type d -name "one_rate_reseeding*")

files=""

for dir in $dirs; do
    name=$(echo $dir | rev | cut -d'/' -f1 | rev)
    files+=" $dir/${name}.log"
done

applauncher NSLogAnalyser -N 1 -noposterior $files -out "$main_dir/one_rate_reseeding_combined.log" > "$main_dir/one_rate_reseeding_combined_terminal.log" 2>&1
