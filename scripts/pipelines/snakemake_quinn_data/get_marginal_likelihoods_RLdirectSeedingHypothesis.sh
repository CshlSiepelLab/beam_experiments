#!/bin/bash

main_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/quinn_2021_lung_cancer_data_2_21_25"


### Combining particles
export applauncher
export main_dir

process_dir() {
    dir=$1
    files=$(find "$dir" -type f -regex '.*/chain_[0-9]+\.log')
    if [[ -n $files ]]; then
        applauncher NSLogAnalyser -N 1 -noposterior $files -out "$dir/combined.log" > "$dir/combined_terminal.log" 2>&1
    fi
}

export -f process_dir

num_threads=50

# gtr
dirs=$(find $main_dir/beam_gtr_ns -maxdepth 2 -mindepth 2 -type d )
echo "$dirs" | parallel -j $num_threads process_dir

# no RL direct seeding
dirs=$(find $main_dir/beam_gtr_noRLdirectSeeding_ns -maxdepth 2 -mindepth 2 -type d )
echo "$dirs" | parallel -j $num_threads process_dir


### Calculating Bayes factors
gtr_dir=$main_dir/beam_gtr_ns

outfile=$main_dir/marginal_likelihoods_noRLdirectSeeding.csv

# Write header
echo "name,gtr_ml,gtr_sd,gtr_noRLdirectSeeding_ml,gtr_noRLdirectSeeding_sd,bf(gtr-gtrNoRLdirectSeeding),diff_threshold" > $outfile

files=$(find $gtr_dir -type f -name "combined_terminal*")

for file in $files; do

    gtr_file=$file
    m2_file=$(echo $gtr_file | sed 's/beam_gtr_ns/beam_gtr_noRLdirectSeeding_ns/g')

    if [[ ! -f $m2_file || ! -f $gtr_file ]]; then
        echo "Missing file: $m2_file or $gtr_file"
        continue
    fi

    name=$(echo $file | rev | cut -d'/' -f2-3 | sed 's/ /_/g' | rev)

    gtr_ml=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f3)
    gtr_sd=$(grep "Marginal likelihood" $gtr_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

    m2_ml=$(grep "Marginal likelihood" $m2_file | cut -d' ' -f3)
    m2_sd=$(grep "Marginal likelihood" $m2_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

    if [[ -z $gtr_ml || -z $m2_ml ]]; then
        beam_bf=""
    else
        beam_bf=$(echo "scale=10; $gtr_ml - $m2_ml" | bc -l)
    fi

    # Required Bayes factor difference threshold based on the estimated standard deviations of the marginal likelihoods from nested sampling
    if [[ -z $gtr_sd || -z $m2_sd ]]; then
        diff_threshold=""
    else
        diff_threshold=$(echo "scale=10; 2 * sqrt(($gtr_sd^2) + ($m2_sd^2))" | bc -l)
    fi
    # Write to outfile
    echo -e "$name,$gtr_ml,$gtr_sd,$m2_ml,$m2_sd,$beam_bf,$diff_threshold" >> $outfile

done

# sort results by Bayes factor from high to low
(head -n1 $outfile &&  tail -n +2 $outfile | sort -t, -k6,6nr)  > $outfile.tmp
mv $outfile.tmp $outfile

# find which need more particles
awk -F',' 'NR > 1 { if (sqrt((5 - $6)^2) < $7) print $1 }' $outfile