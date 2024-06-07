#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime
pipeline_run_name="asymTissueModel_200cell_sims_joint_inference_cassiopeia_machina_5_29_24"
mkdir ${pipeline_run_name}

# copy machina datasets to working directory for the run
dir_pre="data/sim_data_200cellTrees_barcodes_5_27_24"
cp -r ${dir_pre}/* ${pipeline_run_name}/

# intialize files to track metrics for the entire run
accuracy_file="${pipeline_run_name}/accuracy.csv"
echo "dir_name,machina_f1,beast_mcc_f1,beast_posterior_f1,random_f1,consensus_f1,beast_95ci_binary_downsampled_mg,true_mg_machina_f1,true_mg_beast_mcc_f1,true_mg_beast_posterior_f1,true_mg_random_f1,true_mg_consensus_f1,beast_95ci_binary_true_mg,ess_convergence,migration_count,cas_rf_dist,joint_rf_dist,shannon_mut_matrix" > $accuracy_file

commands=()

for dir in ${pipeline_run_name}/*/*;
do
cmd="scripts/pipelines/joint_beast_inference_vs_cassiopeia_machina.sh $dir $accuracy_file"
commands+=("$cmd")
done

for command in "${commands[@]}"
do
  echo "${command}" >> "${pipeline_run_name}/parallel.txt"
done

threads=4
parallel -j 20 < "${pipeline_run_name}/parallel.txt"
rm "${pipeline_run_name}/parallel.txt"
