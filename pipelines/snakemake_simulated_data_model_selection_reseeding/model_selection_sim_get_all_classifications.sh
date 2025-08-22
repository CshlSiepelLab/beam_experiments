#!/bin/bash

main_dir="/grid/siepel/home/staklins/stored_results/beam/latest_results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24"

primaryTissue="P"

all_reseeding="/grid/siepel/home/staklins/stored_results/beam/latest_results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24/reseeding_sims.txt"
chosen_no_reseeding="/grid/siepel/home/staklins/stored_results/beam/latest_results/model_selection_reseeding_no_reseeding_5_8_25_variable_rates_data_8_19_24/no_reseeding_sims_randomly_chosen.txt"

### Combining particles
export applauncher
export main_dir

process_dir() {
    dir=$1

    # Check if combined log file already exists
    if [ -f "$dir/combined_terminal.log" ]; then
        return
    fi

    # To clean improper nested sampling runs that did not fully finish
    files=$(find $dir -type f -name "terminal_*")
    for file in $files; do
        if ! grep -q "Done!" "$file"; then
            run_number=$(basename $file | sed 's/terminal_//g' | sed 's/\.log//g')
            dir=$(dirname $file)
            rm "$file"
            rm $dir/${run_number}.*
            rm $dir/chain_${run_number}.*
        fi
    done

    # Combine logs from all chains in the directory
    files=$(find "$dir" -type f -regex '.*/chain_[0-9]+\.log')
    
    if [ -n "$files" ]; then
        echo $dir
        applauncher NSLogAnalyser -N 1 -noposterior $files -out "$dir/combined.log" > "$dir/combined_terminal.log" 2>&1
    fi
}

export -f process_dir

num_threads=100

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




# ### Separate from above
# # Reset come cps
# cps=(
# mig5_mut0025_1251
# mig5_mut005_21812
# mig5_mut005_23638
# mig6_mut0025_12143
# mig6_mut005_30294
# mig6_mut01_23364
# mig6_mut001_22566
# mig5_mut01_6758
# mig5_mut001_7848
# mig5_mut005_2289
# mig5_mut005_28953
# mig5_mut001_974
# mig5_mut01_2534
# mig6_mut01_26589
# mig6_mut001_22110
# mig6_mut01_28748
# mig6_mut005_28340
# mig6_mut01_13401
# mig5_mut01_6232
# mig6_mut0025_25796
# mig5_mut005_6169
# mig6_mut01_32163
# mig5_mut0025_10050
# mig6_mut001_13259
# mig6_mut001_18037
# mig5_mut0025_6233
# mig5_mut0025_2468
# mig5_mut005_1707
# mig6_mut005_9570
# mig5_mut0025_32617
# mig5_mut01_28920
# mig6_mut005_21643
# mig5_mut001_1688
# mig6_mut001_30120
# mig5_mut005_24716
# mig5_mut005_10364
# mig6_mut0025_1342
# mig6_mut001_10537
# mig5_mut005_10837
# mig6_mut0025_17715
# mig6_mut01_4374
# mig6_mut001_25214
# mig5_mut0025_19909
# mig5_mut01_6228
# mig6_mut005_1596
# mig6_mut01_21493
# )

# for cp in "${cps[@]}"; do
#     rm ${main_dir}/beam_reseeding_ns/${cp}/combined_terminal.log
#     rm ${main_dir}/beam_reseeding_ns/${cp}/combined.log

#     rm ${main_dir}/beam_no_reseeding_ns/${cp}/combined_terminal.log
#     rm ${main_dir}/beam_no_reseeding_ns/${cp}/combined.log
# done


