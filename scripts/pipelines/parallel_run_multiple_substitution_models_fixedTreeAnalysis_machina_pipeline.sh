#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
# source ~/anaconda3/etc/profile.d/conda.sh

# Necessary line to access conda commands on Evolgen lab server (need to make these the same long term)
source ~/miniconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime
pipeline_run_name="multiple_models_beast_machina_performance_marginal_likelihood_3_14_24"
mkdir ${pipeline_run_name}

# copy machina datasets to working directory for the run
data=(m5 m8)
for dataset in ${data[@]}
do
dir_pre="machina_data/sims/"
dir_name="machina_${dataset}_sim_data"
cp -r ${dir_pre}${dir_name} ${pipeline_run_name}/
cp_dir="${pipeline_run_name}/${dir_name}"
done

# intialize files to track metrics for the entire run
accuracy_file="${pipeline_run_name}/accuracy.tsv"
marginal_likelihood_file="${pipeline_run_name}/marginal_likelihoods.tsv"

commands=()

for dir in ${pipeline_run_name}/*/*;
do
cmd="scripts/pipelines/run_multiple_substitution_models_fixedTreeAnalysis_machina_pipeline.sh $dir $pipeline_run_name $accuracy_file $marginal_likelihood_file"
commands+=("$cmd")
done

for command in "${commands[@]}"
do
  echo "${command}" >> "${pipeline_run_name}/parallel.txt"
done

parallel --progress -j 35 < "${pipeline_run_name}/parallel.txt"
rm "${pipeline_run_name}/parallel.txt"
