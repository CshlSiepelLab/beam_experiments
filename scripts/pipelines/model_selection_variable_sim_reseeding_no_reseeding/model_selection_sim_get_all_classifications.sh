#!/bin/bash

main_dir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24"

primaryTissue="P"

all_reseeding="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24/reseeding_sims.txt"
chosen_no_reseeding="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24/no_reseeding_sims_randomly_chosen.txt"

### Combining particles
export applauncher
export main_dir

process_dir() {
    dir=$1

    # Check if combined log file already exists
    if [ -f "$dir/combined_terminal.log" ]; then
        return
    fi

    files=$(find "$dir" -type f -regex '.*/chain_[0-9]+\.log')
    if [[ -n $files ]]; then
        echo $dir
        applauncher NSLogAnalyser -N 1 -noposterior $files -out "$dir/combined.log" > "$dir/combined_terminal.log" 2>&1
    fi
}

export -f process_dir

num_threads=25

# process all directories to combine particles across all chains
dirs=$(find $main_dir/beam_no_reseeding_ns $main_dir/beam_reseeding_ns -maxdepth 1 -mindepth 1 -type d )
echo "$dirs" | parallel -j $num_threads process_dir


# get bayes factors for all simulations
output_csv="$main_dir/beam_classification_results.csv"
echo "sim_name,true,ml_reseeding,sd_reseeding,ml_no_reseeding,sd_no_reseeding,bf,diff_threshold" > $output_csv

true_data_dir="$main_dir/true_reseeding_classification"

sim_names=$(cat $all_reseeding $chosen_no_reseeding | sort | uniq)

for sim_name in $sim_names; do

    true=$(cat $true_data_dir/${sim_name}/reseeding.txt)
    no_reseeding_file="$main_dir/beam_no_reseeding_ns/$sim_name/combined_terminal.log"
    reseeding_file="$main_dir/beam_reseeding_ns/$sim_name/combined_terminal.log"

    if [ -f $no_reseeding_file ] && [ -f $reseeding_file ]; then

        no_reseeding_ml=$(grep "Marginal likelihood" $no_reseeding_file | cut -d' ' -f3)
        no_reseeding_sd=$(grep "Marginal likelihood" $no_reseeding_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        reseeding_ml=$(grep "Marginal likelihood" $reseeding_file | cut -d' ' -f3)
        reseeding_sd=$(grep "Marginal likelihood" $reseeding_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        # Bayes factor is reported with Hnull as the no reseeding and Halt as one rate reseeding, so a positive Bayes factor value supports reseeding and negative supports no reseeding
        beam_bf=$(echo "scale=10; $reseeding_ml - $no_reseeding_ml" | bc)

        # Required Bayes factor difference threshold based on the estimated standard deviations of the marginal likelihoods from nested sampling
        diff_threshold=$(echo "scale=10; 2 * sqrt(($reseeding_sd^2) + ($no_reseeding_sd^2))" | bc -l)
    else
        beam_bf="nan"
        diff_threshold="nan"
        information="nan"
    fi

    echo "$sim_name,$true,$reseeding_ml,$reseeding_sd,$no_reseeding_ml,$no_reseeding_sd,$beam_bf,$diff_threshold" >> $output_csv
done

# sort file by bf
head -n 1 $output_csv > $output_csv.sorted
tail -n +2 $output_csv | sort -t, -k 7,7nr >> $output_csv.sorted
mv $output_csv.sorted $output_csv
