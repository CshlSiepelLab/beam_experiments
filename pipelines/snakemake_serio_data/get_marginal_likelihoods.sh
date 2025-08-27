#!/bin/bash

main_dir="/grid/siepel/home/staklins/stored_results/beam/latest_results/serio_prostate_cancer_data"


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

num_threads=32

# process all chains in one
dirs=$(find $main_dir/beam_random_ns $main_dir/beam_gtr_ns -maxdepth 2 -mindepth 2)
# dirs=$(find $main_dir/beam_reseeding_ns $main_dir/beam_no_reseeding_ns -maxdepth 2 -mindepth 2)

echo "$dirs" | parallel -j $num_threads process_dir



### Calculating Bayes factors
gtr_dir=$main_dir/beam_gtr_ns

outfile=$main_dir/marginal_likelihoods_gtr_random.csv
# outfile=$main_dir/marginal_likelihoods_reseeding_no_reseeding.csv

echo "name,gtr_ml,gtr_sd,random_ml,random_sd,bf(gtr-random),diff_threshold" > $outfile
# echo "name,reseeding_ml,reseeding_sd,no_reseeding_ml,no_reseeding_sd,bf(reseeding-no_reseeding),diff_threshold_reseeding" > $outfile

# # For all cps
files=$(find $gtr_dir -type f -name "combined_terminal*")

# # For select cps
# cps=(
# MMUS1457/CP01
# MMUS1457/CP02
# )

# fileDirs=""
# for cp in "${cps[@]}"; do
#     fileDirs+=" $gtr_dir/$cp"
# done
# files=$(find $fileDirs -type f -name "combined_terminal*")


export outfile

process_file() {
    file=$1

    gtr_file=$file
    random_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_random_ns/g')
    reseeding_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_reseeding_ns/g')
    no_reseeding_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_no_reseeding_ns/g')

    name=$(echo $file | rev | cut -d'/' -f2-3 | sed 's/ /_/g' | rev)

    # Initialize all variables as empty
    gtr_ml=""; gtr_sd=""; random_ml=""; random_sd=""; beam_bf=""; diff_threshold=""
    reseeding_ml=""; reseeding_sd=""; no_reseeding_ml=""; no_reseeding_sd=""; reseeding_bf=""; diff_threshold_reseeding=""

    # GTR vs Random comparison
    if [[ -f $gtr_file && -f $random_file ]]; then
        gtr_ml=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f3)
        gtr_sd=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        random_ml=$(grep "Marginal likelihood" $random_file | cut -d' ' -f3)
        random_sd=$(grep "Marginal likelihood" $random_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        if [[ -n $gtr_ml && -n $random_ml ]]; then
            beam_bf=$(echo "scale=10; $gtr_ml - $random_ml" | bc -l)
        fi

        if [[ -n $gtr_sd && -n $random_sd ]]; then
            diff_threshold=$(echo "scale=10; 2 * sqrt(($gtr_sd^2) + ($random_sd^2))" | bc -l)
        fi
    fi

    # # Reseeding vs No Reseeding comparison
    # if [[ -f $reseeding_file && -f $no_reseeding_file ]]; then
    #     reseeding_ml=$(grep "Marginal likelihood" $reseeding_file | cut -d' ' -f3)
    #     reseeding_sd=$(grep "Marginal likelihood" $reseeding_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

    #     no_reseeding_ml=$(grep "Marginal likelihood" $no_reseeding_file | cut -d' ' -f3)
    #     no_reseeding_sd=$(grep "Marginal likelihood" $no_reseeding_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

    #     if [[ -n $reseeding_ml && -n $no_reseeding_ml ]]; then
    #         reseeding_bf=$(echo "scale=10; $reseeding_ml - $no_reseeding_ml" | bc -l)
    #     fi

    #     if [[ -n $reseeding_sd && -n $no_reseeding_sd ]]; then
    #         diff_threshold_reseeding=$(echo "scale=10; 2 * sqrt(($reseeding_sd^2) + ($no_reseeding_sd^2))" | bc -l)
    #     fi
    # fi

    while ! echo -e "$name,$gtr_ml,$gtr_sd,$random_ml,$random_sd,$beam_bf,$diff_threshold" >> $outfile; do
        sleep 1
    done

    # while ! echo -e "$name,$reseeding_ml,$reseeding_sd,$no_reseeding_ml,$no_reseeding_sd,$reseeding_bf,$diff_threshold_reseeding" >> $outfile; do
    #     sleep 1
    # done
}

export -f process_file

echo "$files" | parallel -j $num_threads process_file


# sort results by Bayes factor from high to low
bf_field=6
(head -n1 $outfile &&  tail -n +2 $outfile | sort -t, -k${bf_field},${bf_field}nr)  > $outfile.tmp
mv $outfile.tmp $outfile






### Separate from above
# To clean combined terminal files for specific mouse_cp combinations

dirs=()

for combo in "${cps[@]}"; do
    mouse=$(echo $combo | cut -d'/' -f1)
    cp=$(echo $combo | cut -d'/' -f2)

    # gtr and random
    gtr_file="$main_dir/beam_gtr_ns/$mouse/$cp/combined_terminal.log"
    random_file="$main_dir/beam_random_ns/$mouse/$cp/combined_terminal.log"

    if [ -f "$gtr_file" ]; then
        rm $gtr_file
        rm $(dirname $gtr_file)/combined.log
    fi
    if [ -f "$random_file" ]; then
        rm "$random_file"
        rm $(dirname $random_file)/combined.log
    fi

    dirs+=("$(dirname $gtr_file)")
    dirs+=("$(dirname $random_file)")

    # # reseeding and no reseeding
    # reseeding_file="$main_dir/beam_reseeding_ns/$mouse/$cp/combined_terminal.log"
    # no_reseeding_file="$main_dir/beam_no_reseeding_ns/$mouse/$cp/combined_terminal.log"

    # if [ -f "$reseeding_file" ]; then
    #     rm $reseeding_file
    #     rm $(dirname $reseeding_file)/combined.log
    # fi
    # if [ -f "$no_reseeding_file" ]; then
    #     rm "$no_reseeding_file"
    #     rm $(dirname $no_reseeding_file)/combined.log
    # fi

    # dirs+=("$(dirname $reseeding_file)")
    # dirs+=("$(dirname $no_reseeding_file)")
done

printf "%s\n" "${dirs[@]}" | parallel -j $num_threads process_dir

