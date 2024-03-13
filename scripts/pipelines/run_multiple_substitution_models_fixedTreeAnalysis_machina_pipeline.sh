#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
# source ~/anaconda3/etc/profile.d/conda.sh

# Necessary line to access conda commands on Evolgen lab server (need to make these the same long term)
source ~/miniconda3/etc/profile.d/conda.sh

# set paths to necessary packages
beast_path=$(which beast)
treeannotator_path=$(which treeannotator)
metastabayes_jar="../metastabayes/metastabayes.jar"

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime
models=(oneRate threeRates sym asym)
pipeline_run_name="multiple_models_beast_machina_performance_marginalikelihood_3_13_24"
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
marginalLikelihood_file="${pipeline_run_name}/marginal_likelihoods.tsv"



for dir in ${pipeline_run_name}/*/*;
do
conda activate simulate
# get m5 or m8 categoruy of machina sim data in which the seed name resides
datatype=$(echo $dir | awk -F'/' '{print $2}')
# get seed name for sim data
dir_prefix=$(echo $dir | awk -F'/' '{print $3}')
# # Format FixedTreeAnalysis input for BEAST2;  Input is simulated tree and tsv of tissue labels; Output is .tree file and .dat file for tissue mapping
sim_tree="${dir}T_${dir_prefix}_unlabeled_true_tree.nwk"
sim_tissues="${dir}T_${dir_prefix}_tissues.tsv"
python ./scripts/format_xml_template_inputs_fixedTreeAnalysis_from_sim.py ${sim_tree} ${sim_tissues}

# Format template xml for each model
seqfile="${dir}T_${dir_prefix}_unlabeled_true_tree_sequences_formatted_for_xml.txt"
taxafile="${dir}T_${dir_prefix}_unlabeled_true_tree_taxonset_formatted_for_xml.txt"
traitfile="${dir}T_${dir_prefix}_unlabeled_true_tree_traitset_formatted_for_xml.txt"
newickfile="${dir}T_${dir_prefix}_unlabeled_true_tree_newick_formatted_for_xml.txt"
primary_tissue="P"
xml_template="inputs/no_bsvss_template_xml_fixedtreeanalysis_machina_sim_universal.xml"
# Longer chain length for asymmetrical setup since convergence is reached later with more parameters
if [ "$model" = "asym" ]; then
    chainlength=10000000
else
    chainlength=1000000
fi
models=(oneRate threeRates sym asym)
commands=()
beast_trees=()
for model in ${models[@]}; do
scripts/format_template_fixedTreeAnalysis_xml_from_sim.sh ${seqfile} ${taxafile} ${traitfile} ${newickfile} ${xml_template} ${primary_tissue} ${chainlength} ${sym_name}
beast_tree="${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}_tissues.tree"
beast_trees+=("$beast_tree")
if [ "$model" = "sym" ] || [ "$model" = "asym" ]; then
    commands+=("${beast_path} -overwrite -working ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}.xml && ${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}_tissues.trees ${beast_tree}")
else
    commands+=("java -jar ${metastabayes_jar} -overwrite -working ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}.xml && ${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}_tissues.trees ${beast_tree}")
fi
done

# run models in parallel in beast
for command in "${commands[@]}"
do
  echo "${command}" >> "${dir}/parallel.txt"
done

parallel -j 4 < "${dir}/parallel.txt"
rm "${dir}/parallel.txt"

# Prep simulated true tree for MACHINA input files
sim_tree_with_tissues="${dir}/T_${dir_prefix}_tissue_labeled_true_tree.nwk"
machina_dir="${dir}/machina"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${sim_tree_with_tissues} ${machina_dir} ${primary_tissue}
conda deactivate

# Run MACHINA
conda activate machina
sed -i '0,/0/s/0/GL/' ${machina_dir}/*.tree     # Prevents MACHINA segmentation fault due to input formatting
./scripts/machina/run_machina.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue ${primary_tissue} --outdir ${machina_dir}
conda deactivate
# Condense MACHINA output into a labeled tree newick format
conda activate simulate
python ./scripts/machina/post_machina_to_tree.py ${sim_tree_with_tissues} ${machina_dir}/T-P-0.labeling ${machina_dir}
conda deactivate
# Remove intermediate MACHINA output files
mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${dir}/
rm -r ${machina_dir}

# Compare results from BEAST2 FixedTreeAnalysis models and MACHINA against the simulated ground truth
conda activate compare_trees
beast_trees_str=$(IFS=','; echo "${beast_trees[*]}")
machina_tree="${dir}/machina_tree_all_tissue_labels.nwk"
python scripts/calculate_internal_node_label_performance.py ${sim_tree_with_tissues} ${beast_trees_str} ${machina_tree} ${dir}
conda deactivate

# Add sub run output to main output files
sim_accuracy_output="${dir}/compare_machina_beast_internal_node_performance.tsv"
if [ ! -s "${accuracy_file}" ]; then
    # If accuracy_file does not exist or is empty then add the first two lines to include the header
    echo -e "$(sed -n '1,2p' ${sim_accuracy_output})" > ${accuracy_file}
else
    # If accuracy_file exists and is not empty
    echo -e "$(sed -n '2p' ${sim_accuracy_output})" >> ${accuracy_file}
fi

done
done
done

