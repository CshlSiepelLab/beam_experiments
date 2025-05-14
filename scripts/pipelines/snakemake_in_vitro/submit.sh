#!/bin/bash

export REPO_PATH=/grid/siepel/home/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/snakemake_in_vitro

snakemake \
--until BeamPerformance \
--latency-wait 30 \
--use-conda \
--use-singularity \
--singularity-args "--bind $HOME/" \
--snakefile $REPO_PATH/Snakefile \
--configfile $REPO_PATH/config/config.yaml \
--printshellcmds \
--keep-going \
--ignore-incomplete \
--cores 10 \
--jobs 10000 \
--cluster-config $REPO_PATH/config/cluster.yaml \
--cluster 'qsub -cwd -pe threads {cluster.cores} -l m_mem_free={cluster.mem} -l h_rt={cluster.runtime} -o {cluster.logout} -e {cluster.logerror}'