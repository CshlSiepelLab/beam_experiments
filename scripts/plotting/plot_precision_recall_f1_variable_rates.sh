#!/bin/bash

indir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_8_19_24_data_from_8_19_24/raw_data"

# get all files
files=$(find $indir -type f -name "*issue_labeled_tree.nwk")

migs=()
muts=()
files_array=()

# get all rate categories
for file in $files; do
    mig=$(echo $file | grep -oP "mig[0-9]+" | grep -oP "[0-9]+")
    mut=$(echo $file | grep -oP "mut[0-9]+" | grep -oP "[0-9]+")
    migs+=($mig)
    muts+=($mut)
    files_array+=($file)
done

# collapse to unique values
migs=($(echo "${migs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
muts=($(echo "${muts[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

mainOutdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_8_19_24_data_from_8_19_24/precision_recall_curve/variable_rates_specific_plots"
mkdir -p $mainOutdir

scripts="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts"

# make output directories and submit jobs for each combination
primary_tissue="P"

for mig in ${migs[@]}; do
    for mut in ${muts[@]}; do
        outdir="${mainOutdir}/mig${mig}_mut${mut}/"
        mkdir -p $outdir

        # subset files to only those with the current pair
        pair="mig${mig}_mut${mut}"
        echo "Processing $pair"
        true_trees=""
        for file in $files; do
            if [[ $file == *$pair* ]]; then
                true_trees+="$file,"
            fi
        done

        # check if true_trees is not empty
        if [[ -n $true_trees ]]; then
            # submit job
            python $scripts/plot_precision_recall.py $true_trees $primary_tissue $outdir
            python $scripts/plot_f1_score.py $outdir/metrics.csv $outdir
        fi
    done
done

