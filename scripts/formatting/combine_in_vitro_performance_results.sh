#!/bin/bash

inDir="/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/results/in_vitro_data_4_24_25"
outfile=$inDir/combined_performance_metrics.csv


files=$(find $inDir -name "*performance.csv")


echo "scenario,cp,method,threshold,true_positives,false_positives,accuracy,precision,f1" > $outfile

for file in $files; do
    if [ -s "$file" ]; then
        cp=$(basename $(dirname $file))
        scenario=$(basename $(dirname $(dirname $(dirname $file))) | grep -o '[0-9]*$')
        method=$(basename $file | cut -d'_' -f1)
        tail -n +2 "$file" | awk -v cp="$cp" -v scenario="$scenario" -v method="$method" '{{print scenario "," cp "," method "," $0}}' >> $outfile
    fi
done