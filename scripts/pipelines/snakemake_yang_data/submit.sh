#!/bin/bash

export REPO_PATH=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/pipelines/snakemake_yang_data

snakemake \
--until runBeam \
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
--cluster-config $REPO_PATH/config/cluster.yaml \
--cluster 'qsub -cwd -pe threads {cluster.cores} -l m_mem_free={cluster.mem} -l h_rt={cluster.runtime} -o {cluster.logout} -e {cluster.logerror}'
# --rerun-incomplete \
