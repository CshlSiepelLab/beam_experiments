#!/bin/bash

export REPO_PATH=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/snakemake_quinn_data

snakemake \
--use-conda \
--use-singularity \
--singularity-args "--bind $HOME/" \
--snakefile $REPO_PATH/Snakefile \
--configfile $REPO_PATH/config/config.yaml \
--printshellcmds \
--keep-going \
--ignore-incomplete \
--cores 1 \
--jobs 500 \
--latency-wait 300 \
--cluster 'qsub -cwd -pe threads {resources.cores} -l m_mem_free={resources.mem_mb} -l h_rt={resources.runtime} -o {log.out} -e {log.err}'
# --rerun-incomplete \
