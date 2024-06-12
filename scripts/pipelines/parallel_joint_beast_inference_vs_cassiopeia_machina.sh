#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime
pipeline_run_name="individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24"
# mkdir ${pipeline_run_name}

# copy machina datasets to working directory for the run
dir_pre="data/new_simulator_unifromTransitionProbs_6_6_24"
# cp -r ${dir_pre}/* ${pipeline_run_name}/

# # NEW ADDITION: send off command to run proper joint inference across datasets
# qsub -cwd -l m_mem_free=5G -pe threads 5 scripts/pipelines/proper_joint_beast_inference.sh $pipeline_run_name

# intialize files to track metrics for the entire run
accuracy_file="${pipeline_run_name}/individual_joint_inference_accuracy.csv"
echo "dir_name,machina_f1,beast_mcc_f1,beast_posterior_f1,random_f1,consensus_f1,beast_95ci_binary_downsampled_mg,true_mg_machina_f1,true_mg_beast_mcc_f1,true_mg_beast_posterior_f1,true_mg_random_f1,true_mg_consensus_f1,beast_95ci_binary_true_mg,ess_convergence,migration_count,cas_rf_dist,joint_rf_dist,shannon_mut_matrix" > $accuracy_file

commands=()
dirs=$(find ${pipeline_run_name}/*/* -maxdepth 0 -type d)
for dir in ${dirs[@]};
do
if [[ $dir == *"proper_joint_beast_inference"* ]]; then
echo "Skip $dir"
else
cmd="scripts/pipelines/joint_beast_inference_vs_cassiopeia_machina.sh $dir $accuracy_file"
commands+=("$cmd")
fi
done

i=0
for command in "${commands[@]}"
do
  # echo "${command}" >> "${pipeline_run_name}/parallel.txt"

  echo "#!/bin/bash" > "${pipeline_run_name}/parallel${i}.sh"
  echo "${command}" >> "${pipeline_run_name}/parallel${i}.sh"
  chmod u+x "${pipeline_run_name}/parallel${i}.sh"
  qsub -cwd -l m_mem_free=5G -pe threads 5 "${pipeline_run_name}/parallel${i}.sh"
  i=$((i+1))
done

# parallel -j 29 < "${pipeline_run_name}/parallel.txt"
# rm "${pipeline_run_name}/parallel.txt"

rm ${pipeline_run_name}/parallel*
