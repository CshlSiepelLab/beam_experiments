#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/model_selection_2_or_3_parameters_reseeding_no_reseeding_12_7_24_variable_rates_data_8_19_24"

primaryTissue="P"

all_reseeding="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/model_selection_full_GTI_reseeding_no_reseeding_12_5_24_variable_rates_data_8_19_24/reseeding_sims.txt"
chosen_no_reseeding="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/model_selection_full_GTI_reseeding_no_reseeding_12_5_24_variable_rates_data_8_19_24/no_reseeding_sims_randomly_chosen.txt"

### MACHINA

data_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24"

outdir_machina="$main_dir/machina_reseeding_classification"
outdir_metient="$main_dir/metient_reseeding_classification"


mkdir -p $outdir_machina

machina_files=$(find $data_dir/machina -type f -name machina_tree_all_tissue_labels.nwk)

for machina_file in $machina_files; do

    sim_id=$(echo $machina_file | rev | cut -d'/' -f2 | rev)

    if ! grep -q $sim_id $all_reseeding && ! grep -q $sim_id $chosen_no_reseeding; then
        continue
    fi

    # evaluate the single graph solution
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/formatting/classify_is_reseeding_from_newick.py $machina_file $primaryTissue $outdir_machina/${sim_id}.txt

done

### METIENT

metient_files=$(find $data_dir/metient -type f -name "*migration_graphs.txt")

mkdir -p $outdir_metient


for metient_file in $metient_files; do

    sim_id=$(echo $metient_file | rev | cut -d'/' -f2 | rev)

    if ! grep -q $sim_id $all_reseeding && ! grep -q $sim_id $chosen_no_reseeding; then
        continue
    fi

    # take the most common classification amongst output graphs
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/formatting/classify_is_reseeding_from_metient.py $metient_file $primaryTissue $outdir_metient/${sim_id}.txt

done

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
echo "sim_name,true,machina,metient,beam_bf,min_bf_diff,information" > $output_csv

true_data_dir="$main_dir/true_reseeding_classification"

sim_names=$(cat $all_reseeding $chosen_no_reseeding | sort | uniq)

for sim_name in $sim_names; do

    if [ -f $true_data_dir/${sim_name}/reseeding.txt ]; then
        true=$(cat $true_data_dir/${sim_name}/reseeding.txt)
    else
        true="nan"
    fi

    if [ -f $outdir_machina/${sim_name}.txt ]; then
        machina=$(cat $outdir_machina/${sim_name}.txt)
    else
        machina="nan"
    fi

    if [ -f $outdir_metient/${sim_name}.txt ]; then
        metient=$(cat $outdir_metient/${sim_name}.txt)
    else
        metient="nan"
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

    echo "$sim_name,$true,$machina,$metient,$beam_bf,$diff_threshold,$information" >> $output_csv
done
