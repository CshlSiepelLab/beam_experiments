#!/bin/bash

# Get all performance threshold data from previous runs

dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates"

outfile="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/bin_information_performance/concat_all_threshold_stats.csv"

files=$(find $dir -type f -name "all_threshold_stats.csv")

echo "Threshold,precision,recall,sim,thresh_counts" > $outfile

for file in $files; do
    tail -n +2 $file >> $outfile
done

# Plot performance related to information

information_csv="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/beam_information_content.csv"

outpdf="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/information_content_GTI_model_12_13_24_variable_rates_data_8_19_24/bin_information_performance/bin_information_precision_recall.pdf"

python ~/bayesian_phylogenetic_metastasis/scripts/plotting/plot_performance_by_information_content_bins.py $outfile $information_csv $outpdf

