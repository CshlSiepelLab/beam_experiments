#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3" 


### Combining particles
export applauncher
export main_dir

process_dir() {
    dir=$1
    files=$(find "$dir" -type f -regex '.*/chain_[0-9]+\.log')
    applauncher NSLogAnalyser -N 1 -noposterior $files -out "$dir/combined.log" > "$dir/combined_terminal.log" 2>&1
}

export -f process_dir

num_threads=40

# gtr
dirs=$(find $main_dir/beam_gtr_ns -maxdepth 2 -mindepth 2 -type d )
echo "$dirs" | parallel -j $num_threads process_dir

# random
dirs=$(find $main_dir/beam_random_ns -maxdepth 2 -mindepth 2 -type d )
echo "$dirs" | parallel -j $num_threads process_dir


### Calculating Bayes factors
gtr_dir=$main_dir/beam_gtr_ns

outfile=$main_dir/marginal_likelihoods.csv

echo "name,gtr_ml,gtr_sd,random_ml,random_sd,bf(gtr-random),diff_threshold" > $outfile

files=$(find $gtr_dir -type f -name "combined_terminal*")

export outfile

process_file() {
    file=$1

    gtr_file=$file
    random_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_random_ns/g')

    if [[ ! -f $random_file || ! -f $gtr_file ]]; then
        echo "Missing file: $random_file or $gtr_file"
        return
    fi

    name=$(echo $file | rev | cut -d'/' -f2-3 | sed 's/ /_/g' | rev)

    gtr_ml=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f3)
    gtr_sd=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

    random_ml=$(grep "Marginal likelihood" $random_file | cut -d' ' -f3)
    random_sd=$(grep "Marginal likelihood" $random_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

    if [[ -z $gtr_ml || -z $random_ml ]]; then
        beam_bf=""
    else
        beam_bf=$(echo "scale=10; $gtr_ml - $random_ml" | bc -l)
    fi

    if [[ -z $gtr_sd || -z $random_sd ]]; then
        diff_threshold=""
    else
        diff_threshold=$(echo "scale=10; 2 * sqrt(($gtr_sd^2) + ($random_sd^2))" | bc -l)
    fi

    while ! echo -e "$name,$gtr_ml,$gtr_sd,$random_ml,$random_sd,$beam_bf,$diff_threshold" >> $outfile; do
        sleep 1
    done
}

export -f process_file

echo "$files" | parallel -j $num_threads process_file


# sort results by Bayes factor from high to low
(head -n1 $outfile &&  tail -n +2 $outfile | sort -t, -k6,6nr)  > $outfile.tmp
mv $outfile.tmp $outfile
