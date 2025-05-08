#!/bin/bash

main_dir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_2_or_3_parameters_reseeding_no_reseeding_12_7_24_variable_rates_data_8_19_24"

primaryTissue="P"

all_reseeding="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_2_or_3_parameters_reseeding_no_reseeding_12_7_24_variable_rates_data_8_19_24/reseeding_sims.txt"
chosen_no_reseeding="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_2_or_3_parameters_reseeding_no_reseeding_12_7_24_variable_rates_data_8_19_24/no_reseeding_sims_randomly_chosen.txt"

data_dir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24"

### BEAM

export applauncher
export main_dir

process_dir() {
    dir=$1
    files=$(find "$dir" -type f -regex '.*/[0-9]+\.log')
    applauncher NSLogAnalyser -N 1 -noposterior $files -out "$dir/combined.log" > "$dir/combined_terminal.log" 2>&1
}

export -f process_dir

num_threads=25

# no rate
dirs=$(find $main_dir/beam_no_reseeding -maxdepth 1 -mindepth 1 -type d )
echo "$dirs" | parallel -j $num_threads process_dir

# one rate
dirs=$(find $main_dir/beam_one_rate_reseeding -maxdepth 1 -mindepth 1 -type d )
echo "$dirs" | parallel -j $num_threads process_dir

### combine classification results to one csv

output_csv="$main_dir/all_classification_results.csv"
echo "sim_name,true,beam_bf,min_bf_diff" > $output_csv

true_data_dir="$main_dir/true_reseeding_classification"

sim_names=$(cat $all_reseeding $chosen_no_reseeding | sort | uniq)

for sim_name in $sim_names; do

    if [ -f $true_data_dir/${sim_name}/reseeding.txt ]; then
        true=$(cat $true_data_dir/${sim_name}/reseeding.txt)
    else
        true="nan"
    fi

    no_results_file="$main_dir/beam_no_reseeding/$sim_name/combined_terminal.log"
    one_results_file="$main_dir/beam_one_rate_reseeding/$sim_name/combined_terminal.log"

    if [ -f $no_results_file ] && [ -f $one_results_file ]; then
        information=$(grep "Marginal likelihood" $no_results_file | cut -d' ' -f6)

        no_ml=$(grep "Marginal likelihood" $no_results_file | cut -d' ' -f3)
        no_sd=$(grep "Marginal likelihood" $no_results_file | cut -d' ' -f4 | cut -d'=' -f5 | sed 's/(//g' | sed 's/)//g')

        one_ml=$(grep "Marginal likelihood" $one_results_file | cut -d' ' -f3)
        one_sd=$(grep "Marginal likelihood" $one_results_file | cut -d' ' -f4 | cut -d'=' -f5 | sed 's/(//g' | sed 's/)//g')

        # Bayes factor is reported with Hnull as the no reseeding and Halt as one rate reseeding, so a positive Bayes factor value supports reseeding and negative supports no reseeding
        beam_bf=$(echo "scale=10; $one_ml - $no_ml" | bc)

        # Required Bayes factor difference threshold based on the estimated standard deviations of the marginal likelihoods from nested sampling
        diff_threshold=$(echo "2 * sqrt(($one_sd^2) + ($no_sd^2))" | bc -l)
    else
        beam_bf="nan"
        diff_threshold="nan"
        information="nan"
    fi

    echo "$sim_name,$true,$beam_bf,$diff_threshold" >> $output_csv
done
