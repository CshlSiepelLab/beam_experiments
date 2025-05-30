#!/bin/bash

for cp in $(seq 11 30); do
    for i in $(seq 1 3); do
        terminal_files=$(find /grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25/beam_gtr${i}/inVitro/CP0${cp} -type f -name "terminal*")
        ess_above_200=false
        for file in $terminal_files; do
            ess_convergence=$(grep "Msamples" $file | tail -n 1 | sed 's/^ *//;s/ *$//' | tr -s ' ' | cut -d' ' -f3)
            if (( $(echo "$ess_convergence > 200" | bc -l) )); then
            ess_above_200=true
        fi
        done
        if [[ "$ess_above_200" == "false" ]]; then
            echo $cp $i
            rm $terminal_files
            
            dir_to_remove=$(dirname $terminal_files)
            find $dir_to_remove -type f -name "*.pdf" -delete
            find $dir_to_remove -type f -name "*.csv" -delete
            find $dir_to_remove -type f -name "*.txt" -delete
            find $dir_to_remove -type f -name "*.pkl" -delete
        fi
    done
done

