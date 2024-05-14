#!/bin/bash

dir="results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/"
files=$(find $dir -type f -name joint_inference_beast_terminal_time.log)

outputfile="results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/convergence_ess.csv"

for file in $files; do
    if [[ $(grep "Operator" $file | wc -l) -ne 0 ]]; then
        ess=$(awk '/Operator/{exit} 1' $file | awk 'NF > 0' | tail -n 1 | awk '{print $3}')
    else
        ess=$(awk 'NF > 0' $file | tail -n 1 | awk '{print $3}')
    fi
    # echo $file
    # echo $ess
    echo $file,$ess >> $outputfile
done