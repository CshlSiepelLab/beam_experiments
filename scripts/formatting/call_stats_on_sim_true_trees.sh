#!/bin/bash

# main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24"
main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24"


primary_tissue="P"
outfile="$main_dir/true_tree_stats.txt"

echo -e "sim_name,migration_count,comigration_count,num_multiedges,met_to_met,reseeding,clonality" > $outfile

# Get all of the true tree labeled newick files to classify the true migration graphs
files=$(find $main_dir/raw_data -type f -name "tissue_labeled_tree.nwk")

for file in $files; do

    sim_name=$(dirname $file | rev | cut -d'/' -f1 | rev)

    # get stats
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/formatting/call_migration_comigration_multiedges_topology_from_labeled_newick.py $sim_name $file $primary_tissue $outfile
done

