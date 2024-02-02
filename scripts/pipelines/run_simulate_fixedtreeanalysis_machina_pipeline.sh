#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
source ~/anaconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime

pipeline_run_name="compare_beast_machina_fixedtree_2_1_24"
mkdir ${pipeline_run_name}

# Specify the number of simulations to be run for ground truth trees with migration data
num_trees=1

for ((i=1; i<=${num_trees}; i++))
do
# Simulate data
conda activate simulate
./scripts/sim_wrapper.sh --design RANDOM --out sim${i} --sites 10 --mutrate 1.0 --samples 50 --migrationrate 1.0 --migration inputs/test_migration_prob_matrix.csv

# Keep only tree and tissue dict, remove barcode data since it is unused for fixed tree analysis
mkdir ${pipeline_run_name}/sim_results_sim${i}
mv sim_results_sim${i}/sim${i}_true.nwk ${pipeline_run_name}/sim_results_sim${i}/
mv sim_results_sim${i}/sim${i}_true_tissues.nwk ${pipeline_run_name}/sim_results_sim${i}/
mv sim_results_sim${i}/sim${i}_tissues.tsv ${pipeline_run_name}/sim_results_sim${i}/
rm -r sim_results_sim${i}

# # Format FixedTreeAnalysis input for BEAST2;  Input is simulated tree and tsv of tissue labels; Output is .tree file and .dat file for tissue mapping
sim_tree="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true.nwk"
sim_tissues="${pipeline_run_name}/sim_results_sim${i}/sim${i}_tissues.tsv"
# python ./scripts/format_fixed_tree_from_sim.py ${sim_tree} ${sim_tissues}
python ./scripts/format_xml_template_inputs_fixedTreeAnalysis_from_sim.py ${sim_tree} ${sim_tissues}

# Format template xml
seqfile="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true_sequences_formatted_for_xml.txt"
taxafile="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true_taxonset_formatted_for_xml.txt"
traitfile="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true_traitset_formatted_for_xml.txt"
newickfile="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true_newick_formatted_for_xml.txt"
scripts/format_template_symmetrical_fixedTreeAnalysis_xml_from_sim.sh ${seqfile} ${taxafile} ${traitfile} ${newickfile}

# Run BEAST2 on formatted xml with output automatically in sim directory
beast_path=$(which beast)
${beast_path} -working ${pipeline_run_name}/sim_results_sim${i}/sim${i}_true_final_input_xml.xml

# Get Maximum Clade Credibility tree from posterior of trees
treeannotator_path=$(which treeannotator)
${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${pipeline_run_name}/sim_results_sim${i}/tissue_tree_with_trait.trees ${pipeline_run_name}/sim_results_sim${i}/tissue_tree_with_trait.tree

# Run MACHINA on simulated true tree
sim_tree_with_tissues="${pipeline_run_name}/sim_results_sim${i}/sim${i}_true_tissues.nwk"
machina_dir="${pipeline_run_name}/sim_results_sim${i}/machina"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${sim_tree_with_tissues} ${machina_dir}
conda deactivate

conda activate machina
./scripts/machina/run_machina.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue t1 --outdir ${machina_dir}
conda deactivate

conda activate simulate
python ./scripts/machina/post_machina_to_tree.py ${tree_dir}/only_leaf_tissue_labels.nwk ${machina_dir}/T-t1-0.labeling ${machina_dir}

# Compare results from BEAST2 FixedTreeAnalysis and MACHINA against the simulated ground truth


conda deactivate
done


