#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
# source ~/anaconda3/etc/profile.d/conda.sh

# Necessary line to access conda commands on Evolgen lab server (need to make these the same long term)
source ~/miniconda3/etc/profile.d/conda.sh

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime


# data=(m5 m8)

# for dataset in ${data[@]}
# do
dataset="m8"

pipeline_run_name="test_fixed_round2_machina_${dataset}_sims_compare_beast_machina_fixedtreeanalysis_default_2_8_24"
mkdir ${pipeline_run_name}

accuracy_file="${pipeline_run_name}/accuracy.tsv"
echo -e "data_id\tmachina\tbeast_strict\tbeast_relaxed\tmachina_nonprimary\tbeast_strict_nonprimary\tbeast_relaxed_nonprimary" > ${accuracy_file}

runtime_file="${pipeline_run_name}/runtime.tsv"
echo -e "data_id\tmachina_seconds\tbeast_seconds" > ${runtime_file}


for dir in machina_${dataset}_sim_data/*/;
do

# dir="machina_m8_sim_data/seed172/"

conda activate simulate
dir_prefix=$(echo $dir | awk -F'/' '{print $2}')


# # Format FixedTreeAnalysis input for BEAST2;  Input is simulated tree and tsv of tissue labels; Output is .tree file and .dat file for tissue mapping
sim_tree="${dir}T_${dir_prefix}_unlabeled_true_tree.nwk"
sim_tissues="${dir}T_${dir_prefix}_tissues.tsv"
# python ./scripts/format_fixed_tree_from_sim.py ${sim_tree} ${sim_tissues}
python ./scripts/format_xml_template_inputs_fixedTreeAnalysis_from_sim.py ${sim_tree} ${sim_tissues}

# Format template xml
seqfile="${dir}T_${dir_prefix}_unlabeled_true_tree_sequences_formatted_for_xml.txt"
taxafile="${dir}T_${dir_prefix}_unlabeled_true_tree_taxonset_formatted_for_xml.txt"
traitfile="${dir}T_${dir_prefix}_unlabeled_true_tree_traitset_formatted_for_xml.txt"
newickfile="${dir}T_${dir_prefix}_unlabeled_true_tree_newick_formatted_for_xml.txt"
primary_tissue="P"
# if [[ $dir == *m5* ]]; then
#     xml_template="inputs/template_xml_symmetrical_machina_sim_m5_data.xml"
# elif [[ $dir == *m8* ]]; then
#     xml_template="inputs/template_xml_symmetrical_machina_sim_m8_data.xml"
# fi
xml_template="inputs/template_xml_fixedtreeanalysis_machina_sim_universal.xml"
symmetric="true"
scripts/format_template_fixedTreeAnalysis_xml_from_sim.sh ${seqfile} ${taxafile} ${traitfile} ${newickfile} ${xml_template} ${primary_tissue} ${symmetric}

# Run BEAST2 on formatted xml with output automatically in sim directory
beast_path=$(which beast)
start_time=$(date +%s.%N)
${beast_path} -overwrite -working ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml.xml
end_time=$(date +%s.%N)
beast_time=$(printf "%.2f" $(echo "$end_time - $start_time" | bc))

# Get Maximum Clade Credibility tree from posterior of trees
treeannotator_path=$(which treeannotator)
${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${dir}/tissue_tree_with_trait.trees ${dir}/tissue_tree_with_trait.tree

# Prep simulated true tree for MACHINA input files
sim_tree_with_tissues="${dir}/T_${dir_prefix}_tissue_labeled_true_tree.nwk"
machina_dir="${dir}/machina"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${sim_tree_with_tissues} ${machina_dir} ${primary_tissue}
conda deactivate

# Run MACHINA
conda activate machina
sed -i '0,/0/s/0/GL/' ${machina_dir}/*.tree     # Prevents MACHINA segmentation fault due to input formatting
start_time=$(date +%s.%N)
./scripts/machina/run_machina.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue ${primary_tissue} --outdir ${machina_dir}
end_time=$(date +%s.%N)
machina_time=$(printf "%.2f" $(echo "$end_time - $start_time" | bc))
conda deactivate

# Condense MACHINA output into a labeled tree newick format
conda activate simulate
python ./scripts/machina/post_machina_to_tree.py ${sim_tree_with_tissues} ${machina_dir}/T-P-0.labeling ${machina_dir}
conda deactivate

# Remove intermediate MACHINA output files
mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${dir}/
rm -r ${machina_dir}

# Compare results from BEAST2 FixedTreeAnalysis and MACHINA against the simulated ground truth
conda activate compare_trees
beast_tree="${dir}/tissue_tree_with_trait.tree"
machina_tree="${dir}/machina_tree_all_tissue_labels.nwk"
python scripts/calculate_internal_node_label_performance.py ${sim_tree_with_tissues} ${beast_tree} ${machina_tree}
conda deactivate

sim_accuracy_output="${dir}/compare_machina_beast_internal_node_performance.tsv"
echo -e "$(sed -n '2p' ${sim_accuracy_output})" >> ${accuracy_file}
echo -e "${dir_prefix}\t${machina_time}\t${beast_time}" >> ${runtime_file}

done
# done

