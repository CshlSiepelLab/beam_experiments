#!/bin/bash

noReseedingTemplate="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/model_selection_mmus1875_cp07/beam_ns_template_no_reseeding_mmus1875_cp07.xml"
oneRateReseedingTemplate="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/model_selection_mmus1875_cp07/beam_ns_template_one_rate_reseeding_mmus1875_cp07.xml"

numReplicates=10

outdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/beam_ns_mmus1875_cp07_no_reseeding_one_rate_reseeding_12_18_24"

mkdir -p $outdir

for i in $(seq 1 $numReplicates); do

    working_dir_no=$outdir/no_reseeding_rep${i}
    mkdir -p $working_dir_no
    cp $noReseedingTemplate $working_dir_no/beam_ns_mmus1875_cp07_no_reseeding_replicate_${i}.xml

    echo -e "java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar -threads 1 -overwrite -working -D outname=no_reseeding_rep${i} $working_dir_no/beam_ns_mmus1875_cp07_no_reseeding_replicate_${i}.xml > $working_dir_no/no_reseeding_rep${i}_terminal.log" >> $outdir/parallel.sh

    working_dir_one=$outdir/one_rate_reseeding_rep${i}
    mkdir -p $working_dir_one
    cp $oneRateReseedingTemplate $working_dir_one/beam_ns_mmus1875_cp07_one_rate_reseeding_replicate_${i}.xml

    echo -e "java -Xmx10g -jar /grid/siepel/home_norepl/staklins/beam/beam.jar -threads 1 -overwrite -working -D outname=one_rate_reseeding_rep${i} $working_dir_one/beam_ns_mmus1875_cp07_one_rate_reseeding_replicate_${i}.xml > $working_dir_one/one_rate_reseeding_rep${i}_terminal.log" >> $outdir/parallel.sh

done

chmod +x $outdir/parallel.sh
parallel --progress -j 20 < $outdir/parallel.sh

