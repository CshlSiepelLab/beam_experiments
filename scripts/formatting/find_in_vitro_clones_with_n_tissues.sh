#!/bin/bash

inDir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25/cp_split"

for cp in $(find $inDir -type f -name "temp_tissues.tsv"); do
    numTissues=$(tail -n +2 $cp | cut -f2 | tr ',' '\n' | grep -v "Pd3" | sort | uniq | wc -l)
    cp_num=$(basename $(dirname $cp))
    if [ $numTissues -eq 7 ]; then
        echo $cp_num $numTissues
    fi
done