#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

# set paths to necessary packages
beast_path=$(which beast)
treeannotator_path=$(which treeannotator)
metastabayes_jar="../metastabayes/metastabayes.jar"

# user input to run from parallel approach
dir=$1
pipeline_run_name=$2
accuracy_file=$3
marginal_likelihood_file=$4

conda activate simulate
# get m5 or m8 categoruy of machina sim data in which the seed name resides
datatype=$(echo $dir | awk -F'/' '{print $2}')
# get seed name for sim data
dir_prefix=$(echo $dir | awk -F'/' '{print $3}')
# # Format FixedTreeAnalysis input for BEAST2;  Input is simulated tree and tsv of tissue labels; Output is .tree file and .dat file for tissue mapping
sim_tree="${dir}/T_${dir_prefix}_unlabeled_true_tree.nwk"
sim_tissues="${dir}/T_${dir_prefix}_tissues.tsv"
python ./scripts/format_xml_template_inputs_fixedTreeAnalysis_from_sim.py ${sim_tree} ${sim_tissues}

# Format template xml for each model
seqfile="${dir}/T_${dir_prefix}_unlabeled_true_tree_sequences_formatted_for_xml.txt"
taxafile="${dir}/T_${dir_prefix}_unlabeled_true_tree_taxonset_formatted_for_xml.txt"
traitfile="${dir}/T_${dir_prefix}_unlabeled_true_tree_traitset_formatted_for_xml.txt"
newickfile="${dir}/T_${dir_prefix}_unlabeled_true_tree_newick_formatted_for_xml.txt"
primary_tissue="P"
xml_template="inputs/no_bsvss_template_xml_fixedtreeanalysis_machina_sim_universal.xml"
beast_trees=()
models=(oneRate threeRates sym asym)
for model in ${models[@]}; do
# Longer chain length for asymmetrical setup since convergence is reached later with more parameters
if [ "$model" = "asym" ]; then
    chainlength=10000000
else
    chainlength=1000000
fi
scripts/format_template_fixedTreeAnalysis_xml_from_sim.sh ${seqfile} ${taxafile} ${traitfile} ${newickfile} ${xml_template} ${primary_tissue} ${chainlength} ${model}
beast_tree="${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}_tissues.tree"
beast_trees+=("$beast_tree")
xml_path="${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}.xml"
ns_dir="${dir}/${model}_nested_sampling"
active_particles=10
subchainlen=10000
if [ "$model" = "sym" ] || [ "$model" = "asym" ]; then
    ${beast_path} -overwrite -working ${xml_path}
    ${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}_tissues.trees ${beast_tree}
    mkdir ${ns_dir}
    scripts/nested_sampling_marginal_likelihood_from_xml.sh --xml ${xml_path} --dir ${ns_dir} --active_particles ${active_particles} --sub_chain_length ${subchainlen}
else
    java -jar ${metastabayes_jar} -overwrite -working ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}.xml
    ${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${dir}/T_${dir_prefix}_unlabeled_true_tree_final_input_xml_${model}_tissues.trees ${beast_tree}
    mkdir ${ns_dir}
    scripts/nested_sampling_marginal_likelihood_from_xml.sh --xml ${xml_path} --dir ${ns_dir} --active_particles ${active_particles} --sub_chain_length ${subchainlen} --model metastabayes
fi
done

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
python scripts/calculate_internal_node_label_performance.py ${sim_tree_with_tissues} ${beast_trees_str} ${machina_tree} ${dir} ${datatype}
conda deactivate

# Add sub run outputs to main output files
sim_accuracy_output="${dir}/compare_machina_beast_internal_node_performance.csv"
if [ ! -s "${accuracy_file}" ]; then
    # If accuracy_file does not exist or is empty then add the first two lines to include the header
    # echo -e "$(sed -n '1,2p' ${sim_accuracy_output})" > ${accuracy_file}
    echo -e "$(sed -n '2p' ${sim_accuracy_output})" >> ${accuracy_file}
else
    # If accuracy_file exists and is not empty
    echo -e "$(sed -n '2p' ${sim_accuracy_output})" >> ${accuracy_file}
fi

if [ ! -s "${marginal_likelihood_file}" ]; then
ml_header="data_id"
for m in ${models[@]}; do
ml_header+=",${m}_ml,${m}_sd"
done
echo -e $ml_header > $marginal_likelihood_file
fi

data_id=$(basename "$dir")
ml_str="${data_id}_${datatype}"

for mod in ${models[@]}; do
ml_output_path="${dir}/${mod}_nested_sampling/xml1/xml1_marginal_likelihood_run.txt"
ml_line=$(grep "Marginal likelihood:" "$ml_output_path" | grep -vE "subsample|bootstrap")
ml=$(echo "$ml_line" | cut -d ' ' -f 3 | cut -d '(' -f 1)
sd=$(echo "$ml_line" | cut -d ' ' -f 3 | cut -d '(' -f 2 | cut -d ')' -f 1)
ml_str+=",${ml},${sd}"
done
echo -e $ml_str >> $marginal_likelihood_file

