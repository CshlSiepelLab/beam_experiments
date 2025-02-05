#!/bin/bash

main_dir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_1_20_25_asv_cutoff_3

# Combining mcmc chains if not already done
logcombiner="/grid/siepel/home_norepl/staklins/bin/beast/bin/logcombiner"
essCutoff=200

export logcombiner essCutoff

combine_logs() {
    dir=$1
    log_files=""
    trees_files=""
    for file in $(find $dir -type f -name "chain_*.log"); do
        terminal_log=$(echo $file | sed 's/chain_/terminal_/')
        ess_convergence=$(grep "Msamples" $terminal_log | tail -n 1 | sed 's/^ *//;s/ *$//' | tr -s ' ' | cut -d' ' -f3)
        if (( $(echo "$ess_convergence > $essCutoff" | bc -l) )); then
            log_files+="-log $file "
            trees_file=$(echo $file | sed 's/\.log$/.trees/')
            trees_files+="-log $trees_file "
        else
            echo "$file did not converge, so not including it in the final output"
        fi
    done
    if [ -n "$log_files" ]; then
        $logcombiner $log_files -o $dir/combined.log
        $logcombiner $trees_files -o $dir/combined.trees
    fi
}

export -f combine_logs

find $main_dir/beam_gtr -mindepth 2 -maxdepth 2 -type d | parallel -j 50 combine_logs {}

# Mutual information part
threads=5

files=$(find $main_dir/beam_gtr -type f -name "combined.trees")
primary_tissue="PRL"

for file in $files; do
    dir=$(dirname $file)
    qsub -cwd -l m_mem_free=1G -pe threads $threads -b y "python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/mutual_information_from_beast_posterior.py $file $primary_tissue $dir $threads"    exit
done

# Collect results
result_files=$(find $main_dir/beam_gtr -type f -name "posterior_trees_migration_mutual_information.txt")

outfile="$main_dir/gtr_beam_mutual_information.csv"

echo "mouse_cp,mutual_information_normalized" > $outfile

for file in $result_files; do
    name=$(dirname $file | rev | cut -d'/' -f1-2 | rev | sed "s/\//_/")
    mutual_information=$(head -n 1 $file | tr -d '[:space:]')
    echo "$name,$mutual_information" >> $outfile
done
