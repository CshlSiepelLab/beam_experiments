#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
source ~/anaconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime

pipeline_run_name="compare_beast_machina_fixedtree_2_1_24"
mkdir ${pipeline_run_name}

conda activate simulate

# Specify the number of simulations to be run for ground truth trees with migration data
num_trees=1

for ((i=1; i<=${num_trees}; i++))
do
# Simulate data
./scripts/sim_wrapper.sh --design RANDOM --out sim${i} --sites 10 --mutrate 1.0 --samples 50 --migration inputs/test_migration_prob_matrix.csv

# Keep only tree and tissue dict, remove barcode data since it is unused for fixed tree analysis
mkdir ${pipeline_run_name}/sim_results_sim${i}
mv sim_results_sim${i}/sim${i}_true.nwk ${pipeline_run_name}/sim_results_sim${i}/
mv sim_results_sim${i}/sim${i}_tissues.tsv ${pipeline_run_name}/sim_results_sim${i}/
rm -r sim_results_sim${i}

# Format FixedTreeAnalysis input for BEAST2;  Input is simulated tree and tsv of tissue labels; Output is .tree file and .dat file for tissue mapping
sim_tree="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true.nwk"
sim_tissues="${pipeline_run_name}/sim_results_sim${i}/sim${i}_tissues.tsv"
python ./scripts/format_fixed_tree_from_sim.py ${sim_tree} ${sim_tissues}

done

conda deactivate
