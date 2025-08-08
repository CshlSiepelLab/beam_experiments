#!/bin/bash

export REPO_PATH=/grid/siepel/home/staklins/projects/crispr_barcode/bayesian_phylogenetic_metastasis/pipelines/snakemake_serio_data

snakemake \
--until runLAML \
--latency-wait 30 \
--use-conda \
--use-singularity \
--singularity-args "--bind $HOME/ --bind $GRB_LICENSE_FILE:/mnt/gurobi.lic --env GRB_LICENSE_FILE=/mnt/gurobi.lic" \
--snakefile $REPO_PATH/Snakefile \
--configfile $REPO_PATH/config/config.yaml \
--printshellcmds \
--keep-going \
--ignore-incomplete \
--cores 1 \
--jobs 10000 \
--cluster-config $REPO_PATH/config/cluster.yaml \
--cluster 'qsub -cwd -pe threads {cluster.cores} -l m_mem_free={cluster.mem} -l h_rt={cluster.runtime} -o {cluster.logout} -e {cluster.logerror}'
