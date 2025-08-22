#!/bin/bash

main_dir="/grid/siepel/home/staklins/stored_results/beam/latest_results/quinn_2021_lung_cancer_data" 


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

# process all chains in one
# dirs=$(find $main_dir/beam_gtr_ns $main_dir/beam_random_ns $main_dir/beam_gtr_noRLdirectSeeding_ns -maxdepth 2 -mindepth 2 -type d)
dirs=$(find $main_dir/beam_gtr_noRLdirectSeeding_ns $main_dir/beam_gtr_ns -maxdepth 2 -mindepth 2 -type d)
echo "$dirs" | parallel -j $num_threads process_dir


### Calculating Bayes factors
gtr_dir=$main_dir/beam_gtr_ns

outfile=$main_dir/marginal_likelihoods.csv

# Write header
echo "name,gtr_ml,gtr_sd,random_ml,random_sd,bf(gtr-random),diff_threshold_random,gtr_noRLdirectSeeding_ml,gtr_noRLdirectSeeding_sd,bf(gtr-gtrNoRLdirectSeeding),diff_threshold_RL" > $outfile

files=$(find $gtr_dir -type f -name "combined_terminal*")

for file in $files; do
    gtr_file=$file
    random_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_random_ns/g')
    noRL_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_gtr_noRLdirectSeeding_ns/g')

    name=$(echo $file | rev | cut -d'/' -f2-3 | sed 's/ /_/g' | rev)
    
    # Initialize all variables as empty
    gtr_ml=""; gtr_sd=""; random_ml=""; random_sd=""; noRL_ml=""; noRL_sd=""
    beam_bf_random=""; diff_threshold_random=""; beam_bf_RL=""; diff_threshold_RL=""

    # First handle GTR and random comparison
    if [[ -f $gtr_file && -f $random_file ]]; then
        gtr_ml=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f3)
        gtr_sd=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        random_ml=$(grep "Marginal likelihood" $random_file | cut -d' ' -f3)
        random_sd=$(grep "Marginal likelihood" $random_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        # Calculate Bayes factor for GTR vs random
        if [[ -n $gtr_ml && -n $random_ml ]]; then
            beam_bf_random=$(echo "scale=10; $gtr_ml - $random_ml" | bc -l)
        fi

        # Calculate threshold for GTR vs random
        if [[ -n $gtr_sd && -n $random_sd ]]; then
            diff_threshold_random=$(echo "scale=10; 2 * sqrt(($gtr_sd^2) + ($random_sd^2))" | bc -l)
        fi
    else
        echo "Missing file: $gtr_file or $random_file"
    fi

    # Then handle RL comparison if available
    if [[ -f $noRL_file ]]; then
        noRL_ml=$(grep "Marginal likelihood" $noRL_file | cut -d' ' -f3)
        noRL_sd=$(grep "Marginal likelihood" $noRL_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        # Calculate Bayes factor for GTR vs noRL
        if [[ -n $gtr_ml && -n $noRL_ml ]]; then
            beam_bf_RL=$(echo "scale=10; $gtr_ml - $noRL_ml" | bc -l)
        fi

        # Calculate threshold for GTR vs noRL
        if [[ -n $gtr_sd && -n $noRL_sd ]]; then
            diff_threshold_RL=$(echo "scale=10; 2 * sqrt(($gtr_sd^2) + ($noRL_sd^2))" | bc -l)
        fi
    fi

    # Write to outfile
    echo -e "$name,$gtr_ml,$gtr_sd,$random_ml,$random_sd,$beam_bf_random,$diff_threshold_random,$noRL_ml,$noRL_sd,$beam_bf_RL,$diff_threshold_RL" >> $outfile
done


# sort results by Bayes factor from high to low
bf_field=6
(head -n1 $outfile &&  tail -n +2 $outfile | sort -t, -k${bf_field},${bf_field}nr)  > $outfile.tmp
mv $outfile.tmp $outfile




### Separate from above
# Reset come cps
cps=($(seq 1 100))

for cp in "${cps[@]}"; do
    if [ $cp -gt 50 ]; then
        rm ${main_dir}/beam_gtr_noRLdirectSeeding_ns/5k/${cp}/combined_terminal.log
        rm ${main_dir}/beam_gtr_noRLdirectSeeding_ns/5k/${cp}/combined.log

        rm ${main_dir}/beam_gtr_ns/5k/${cp}/combined_terminal.log
        rm ${main_dir}/beam_gtr_ns/5k/${cp}/combined.log
    fi
done


