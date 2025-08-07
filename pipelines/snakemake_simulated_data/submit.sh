#!/bin/bash

export REPO_PATH=/grid/siepel/home/staklins/projects/crispr_barcode/bayesian_phylogenetic_metastasis/pipelines/snakemake_simulated_data

snakemake \
--until runLAML \
--use-singularity \
--singularity-args "--bind $HOME/" \
--latency-wait 15 \
--use-conda \
--snakefile $REPO_PATH/Snakefile \
--configfile $REPO_PATH/config/config.yaml \
--printshellcmds \
--keep-going \
--ignore-incomplete \
--cores 1 \
--jobs 100 \
--cluster-config $REPO_PATH/config/cluster.yaml \
--cluster 'qsub -cwd -pe threads {cluster.cores} -l m_mem_free={cluster.mem} -l h_rt={cluster.runtime} -o {cluster.logout} -e {cluster.logerror}'
