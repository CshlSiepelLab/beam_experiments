#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
# source ~/anaconda3/etc/profile.d/conda.sh

# Necessary line to access conda commands on lab server (need to make these the same long term)
source ~/miniconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime

pipeline_run_name="machina_m5_sims_compare_beast_machina_fixedtreeanalysis_default_2_7_24"
mkdir ${pipeline_run_name}

accuracy_file="${pipeline_run_name}/accuracy.tsv"
echo -e "data_id\tmachina\tbeast_strict\tbeast_relaxed\tmachina_nonprimary\tbeast_strict_nonprimary\tbeast_relaxed_nonprimary" > ${accuracy_file}

runtime_file="${pipeline_run_name}/runtime.tsv"
echo -e "data_id\tmachina_seconds\tbeast_seconds" > ${runtime_file}

for dir in machina_m5_sim_data/*/;
do
conda activate simulate
dir_prefix=$(echo $dir | awk -F'/' '{print $2}')
# # Format FixedTreeAnalysis input for BEAST2;  Input is simulated tree and tsv of tissue labels; Output is .tree file and .dat file for tissue mapping
sim_tree="${dir}T_${dir_prefix}_unlabeled_true_tree.nwk"
sim_tissues="${dir}T_${dir_prefix}_tissues_tsv"
# python ./scripts/format_fixed_tree_from_sim.py ${sim_tree} ${sim_tissues}
python ./scripts/format_xml_template_inputs_fixedTreeAnalysis_from_sim.py ${sim_tree} ${sim_tissues}

# Format template xml
seqfile="${dir}T_${dir_prefix}_unlabeled_true_tree_sequences_formatted_for_xml.txt"
taxafile="${dir}T_${dir_prefix}_unlabeled_true_tree_taxonset_formatted_for_xml.txt"
traitfile="${dir}T_${dir_prefix}_unlabeled_true_tree_traitset_formatted_for_xml.txt"
newickfile="${dir}T_${dir_prefix}_unlabeled_true_tree_newick_formatted_for_xml.txt"
xml_template="inputs/template_xml_symmetrical_machina_sim_data.xml"
scripts/format_template_symmetrical_fixedTreeAnalysis_xml_from_sim.sh ${seqfile} ${taxafile} ${traitfile} ${newickfile} ${xml_template}

# Run BEAST2 on formatted xml with output automatically in sim directory
beast_path=$(which beast)
start_time=$(date +%s.%N)
${beast_path} -working ${dir_prefix}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml.xml
end_time=$(date +%s.%N)
beast_time=$(printf "%.2f" $(echo "$end_time - $start_time" | bc))

# Get Maximum Clade Credibility tree from posterior of trees
treeannotator_path=$(which treeannotator)
${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${pipeline_run_name}/${dir_prefix}/tissue_tree_with_trait.trees ${pipeline_run_name}/${dir_prefix}/tissue_tree_with_trait.tree

# Prep simulated true tree for MACHINA input files
sim_tree_with_tissues="${dir_prefix}/T_${dir_prefix}_tissue_labeled_true_tree.nwk"
machina_dir="${pipeline_run_name}/${dir_prefix}/machina"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${sim_tree_with_tissues} ${machina_dir}
conda deactivate

# Run MACHINA
conda activate machina
start_time=$(date +%s.%N)
./scripts/machina/run_machina.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue t1 --outdir ${machina_dir}
end_time=$(date +%s.%N)
machina_time=$(printf "%.2f" $(echo "$end_time - $start_time" | bc))
conda deactivate

# Condense MACHINA output into a labeled tree newick format
conda activate simulate
python ./scripts/machina/post_machina_to_tree.py ${sim_tree_with_tissues} ${machina_dir}/T-t1-0.labeling ${machina_dir}
conda deactivate

# Remove intermediate MACHINA output files
mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${pipeline_run_name}/${dir_prefix}/
rm -r ${machina_dir}

# Compare results from BEAST2 FixedTreeAnalysis and MACHINA against the simulated ground truth
conda activate compare_trees
beast_tree="${pipeline_run_name}/${dir_prefix}/tissue_tree_with_trait.tree"
machina_tree="${pipeline_run_name}/${dir_prefix}/machina_tree_all_tissue_labels.nwk"
python scripts/calculate_internal_node_label_performance.py ${sim_tree_with_tissues} ${beast_tree} ${machina_tree}
conda deactivate

sim_accuracy_output="${pipeline_run_name}/${dir_prefix}/compare_machina_beast_internal_node_performance.tsv"
echo -e "$(sed -n '2p' ${sim_accuracy_output})" >> ${accuracy_file}
echo -e "sim${i}\t${machina_time}\t${beast_time}" >> ${runtime_file}

done


