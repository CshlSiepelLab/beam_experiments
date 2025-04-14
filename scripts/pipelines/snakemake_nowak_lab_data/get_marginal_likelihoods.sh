#!/bin/bash

main_dir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_2_24_25"


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
    
    if [ -n "$files" ]; then
        applauncher NSLogAnalyser -N 1 -noposterior $files -out "$dir/combined.log" > "$dir/combined_terminal.log" 2>&1
    fi
}

export -f process_dir

num_threads=50

# process all chains in one
dirs=$(find $main_dir/beam_gtr_ns $main_dir/beam_random_ns $main_dir/beam_reseeding_ns $main_dir/beam_no_reseeding_ns -maxdepth 2 -mindepth 2)
echo "$dirs" | parallel -j $num_threads process_dir


### Calculating Bayes factors
gtr_dir=$main_dir/beam_gtr_ns

outfile=$main_dir/marginal_likelihoods.csv

echo "name,gtr_ml,gtr_sd,random_ml,random_sd,bf(gtr-random),diff_threshold,reseeding_ml,reseeding_sd,no_reseeding_ml,no_reseeding_sd,bf(reseeding-no_reseeding),diff_threshold_reseeding" > $outfile

files=$(find $gtr_dir -type f -name "combined_terminal*")

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

    # Reseeding vs No Reseeding comparison
    if [[ -f $reseeding_file && -f $no_reseeding_file ]]; then
        reseeding_ml=$(grep "Marginal likelihood" $reseeding_file | cut -d' ' -f3)
        reseeding_sd=$(grep "Marginal likelihood" $reseeding_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        no_reseeding_ml=$(grep "Marginal likelihood" $no_reseeding_file | cut -d' ' -f3)
        no_reseeding_sd=$(grep "Marginal likelihood" $no_reseeding_file | cut -d' ' -f4 | cut -d'(' -f4 | sed 's/)//g')

        if [[ -n $reseeding_ml && -n $no_reseeding_ml ]]; then
            reseeding_bf=$(echo "scale=10; $reseeding_ml - $no_reseeding_ml" | bc -l)
        fi

        if [[ -n $reseeding_sd && -n $no_reseeding_sd ]]; then
            diff_threshold_reseeding=$(echo "scale=10; 2 * sqrt(($reseeding_sd^2) + ($no_reseeding_sd^2))" | bc -l)
        fi
    fi

    while ! echo -e "$name,$gtr_ml,$gtr_sd,$random_ml,$random_sd,$beam_bf,$diff_threshold,$reseeding_ml,$reseeding_sd,$no_reseeding_ml,$no_reseeding_sd,$reseeding_bf,$diff_threshold_reseeding" >> $outfile; do
        sleep 1
    done
}

export -f process_file

echo "$files" | parallel -j $num_threads process_file


# sort results by Bayes factor from high to low
bf_field=6
(head -n1 $outfile &&  tail -n +2 $outfile | sort -t, -k${bf_field},${bf_field}nr)  > $outfile.tmp
mv $outfile.tmp $outfile

# find which need more particles
diff_field=$(( bf_field + 1 ))
threshold=0
awk -F',' -v bf=$bf_field -v diff=$diff_field 'NR > 1 { if (sqrt(($threshold - $bf)^2) < $diff) print $1 }' $outfile

