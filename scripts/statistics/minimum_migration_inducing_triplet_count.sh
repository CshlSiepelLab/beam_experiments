#! /bin/bash

main_dir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25"

primary_tissue="P"

num_scenarios=3

outfile="$main_dir/minimum_migration_inducing_triplet_counts.csv"

echo "scenario,cp,mmitc" > $outfile

for scenario in $(seq 1 $num_scenarios); do
    echo "Scenario $scenario"

    matrix_files=$(find $main_dir/beam_prep$scenario/ -type f -name "expanded_clones_matrix.tsv")

    for matrix_file in $matrix_files; do

        cp="$(basename $(dirname $matrix_file))"
        echo "$cp"

        char_matrix_fp="$matrix_file"
        tissue_labels_fp="$(dirname $matrix_file)/expanded_clones_tissues.tsv"

        count="nan"

        if [ -s $char_matrix_fp ] && [ -s $tissue_labels_fp ]; then
            count=$(python /grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/minimum_migration_inducing_triplet_count.py \
            $char_matrix_fp $tissue_labels_fp $primary_tissue | grep  "Minimum Migration-Inducing Triplet Count:" | cut -d' ' -f 5)
        fi

        echo "$scenario,$cp,$count" >> $outfile
    
    done
done