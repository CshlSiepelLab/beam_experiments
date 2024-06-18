#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime
pipeline_run_name="barcodeSites50_uniform_rates_precision_recall_joint_inference_vs_cassiopeia_machina_6_18_24"
mkdir ${pipeline_run_name}

# copy machina datasets to working directory for the run
dir_pre="data/targetSites50_uniformMigration_6_18_24"
cp -r ${dir_pre}/* ${pipeline_run_name}/

# # NEW ADDITION: send off command to run proper joint inference across datasets
# qsub -cwd -l m_mem_free=5G -pe threads 5 scripts/pipelines/proper_joint_beast_inference.sh $pipeline_run_name

# intialize files to track metrics for the entire run
accuracy_file="${pipeline_run_name}/individual_joint_inference_accuracy.csv"
echo "dir,\
      tree_machina_f1,\
      tree_machina_precision,\
      tree_machina_recall,\
      tree_beast_mcc_f1,\
      tree_beast_mcc_precision,\
      tree_beast_mcc_recall,\
      tree_beast_posterior_f1,\
      tree_beast_posterior_precision,\
      tree_beast_posterior_recall,\
      tree_random_f1,\
      tree_random_precision,\
      tree_random_recall,\
      tree_consensus_f1,\
      tree_consensus_precision,\
      tree_consensus_recall,\
      beast_posterior_95ci_binary,\
      true_machina_f1,\
      true_machina_precision,\
      true_machina_recall,\
      true_beast_mcc_f1,\
      true_beast_mcc_precision,\
      true_beast_mcc_recall,\
      true_beast_posterior_f1,\
      true_beast_posterior_precision,\
      true_beast_posterior_recall,\
      true_random_f1,\
      true_random_precision,\
      true_random_recall,\
      true_consensus_f1,\
      true_consensus_precision,\
      true_consensus_recall,\
      true_beast_posterior_95ci_binary,\
      ess_convergences[@],\
      migration_count,\
      cas_rf_dist,\
      joint_rf_dist,\
      shannon_mut_matrix" > ${accuracy_file}

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

job_logs_dir="${pipeline_run_name}/job_logs"
mkdir $job_logs_dir

i=0
for command in "${commands[@]}"
do
  # echo "${command}" >> "${pipeline_run_name}/parallel.txt"
  cmd_file="${job_logs_dir}/parallel${i}.sh"
  echo "#!/bin/bash" > $cmd_file
  echo "${command}" >> $cmd_file
  chmod u+x $cmd_file
  qsub -cwd -l m_mem_free=5G -pe threads 5  -o ${cmd_file}.out -e ${cmd_file}.err $cmd_file
  i=$((i+1))
done

# parallel -j 29 < "${pipeline_run_name}/parallel.txt"
# rm "${pipeline_run_name}/parallel.txt"

