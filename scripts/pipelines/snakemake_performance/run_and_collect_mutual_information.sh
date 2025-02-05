#!/bin/bash

main_dir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_1_31_25_data_from_8_19_24
threads=5

files=$(find $main_dir/beam_gtr -type f -name "combined.trees")
primary_tissue="P"

for file in $files; do
    dir=$(dirname $file)
    qsub -cwd -l m_mem_free=1G -pe threads $threads -b y "python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/mutual_information_from_beast_posterior.py $file $primary_tissue $dir $threads"    exit
done

# Collect results
result_files=$(find $main_dir/beam_gtr -type f -name "posterior_trees_migration_mutual_information.txt")

outfile="$main_dir/gtr_beam_mutual_information.csv"

echo "simname,mutual_information_normalized" > $outfile

for file in $result_files; do
    simname=$(basename $(dirname $file))
    mutual_information=$(head -n 1 $file | tr -d '[:space:]')
    echo "$simname,$mutual_information" >> $outfile
done
