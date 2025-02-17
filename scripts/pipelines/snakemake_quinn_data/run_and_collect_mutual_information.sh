#!/bin/bash

main_dir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_1_22_25

# Mutual information part
threads=5

files=$(find $main_dir/beam_gtr -type f -name "combined.trees")
primary_tissue="LL"

i=0
for file in $files; do
    dir=$(dirname $file)

    # skip if already computed
    if [ -f $dir/posterior_trees_migration_mutual_information.txt ]; then
        continue
    fi
    
    qsub -cwd -l m_mem_free=1G -pe threads $threads -N "run$i" -b y "python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/mutual_information_from_beast_posterior.py $file $primary_tissue $dir $threads"
    i=$((i+1))
done

# Collect results
result_files=$(find $main_dir/beam_gtr -type f -name "posterior_trees_migration_mutual_information.txt")

outfile="$main_dir/gtr_beam_mutual_information.csv"

echo "mouse_cp,mutual_information_normalized" > $outfile

for file in $result_files; do
    name=$(dirname $file | rev | cut -d'/' -f1-2 | rev | sed "s/\//_/")
    mutual_information=$(head -n 1 $file | tr -d '[:space:]')
    echo "$name,$mutual_information" >> $outfile
done
