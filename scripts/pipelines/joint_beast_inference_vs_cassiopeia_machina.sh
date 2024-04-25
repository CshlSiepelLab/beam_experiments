#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
# source ~/anaconda3/etc/profile.d/conda.sh

### This pipeline takes in simulated data in the form of an indel character matrix and ground truth tree with tissue labels and the compares cassiopeia->machina and joint tree and tissue BEAST inference method for performance in inferring the migration graph vs the ground truth

# # user inputs
directory=$1
accuracy_file=$2

# directory="./11031"
# accuracy_file="./accuracy.csv"

sim_matrix=${directory}/*_indel_character_matrix.tsv
true_tree=${directory}/cell_tree_seed*.nwk
true_tissues=${directory}/cell_tree_seed*.vertex.labeling
leaf_tissues=$(ls ${directory}/cell_tree_*[0-9].labeling)
drivers=${directory}/drivers_seed*.txt


# for testing
# sim_matrix="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/24874_indel_character_matrix.tsv"
# true_tree="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874.nwk"
# true_tissues="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874.vertex.labeling"
# leaf_tissues="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874.labeling"

# get executable paths for beast
treeannotator_path=$(which treeannotator)
metastabayes_jar="../metastabayes/metastabayes.jar"

# get working dir
dir=$(dirname "$sim_matrix")

# run cassiopeia-greedy on matrix
conda activate simulate
python scripts/cassiopeia_greedy.py $sim_matrix

# run machina on cassiopeia-greedy inferred tree
# Prep cas tree for MACHINA input files
cas_tree_tissues="${dir}/cassiopeia_greedy_inferred.nwk"
machina_dir="${dir}/machina"
primary_tissue="P"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${cas_tree_tissues} ${machina_dir} ${primary_tissue} ${leaf_tissues}
conda deactivate

# Run MACHINA
conda activate machina
#sed -i '0,/0/s/0/GL/' ${machina_dir}/*.tree     # Prevents MACHINA segmentation fault due to input formatting
./scripts/machina/run_machina_tr.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue ${primary_tissue} --outdir ${machina_dir}
conda deactivate

# Condense MACHINA output into a labeled tree newick format
conda activate simulate
python ./scripts/machina/post_machina_tr_to_tree.py ${machina_dir}/P-T-P-R.tree ${machina_dir}/P-T-P-R.labeling ${machina_dir}
conda deactivate
# Remove intermediate MACHINA output files
machina_tree="${dir}/machina_tree_all_tissue_labels.nwk" 
mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${machina_tree}
#rm -r ${machina_dir}

# setup joint beast inference
template_xml="inputs/joint_inference_beast_template.xml"
sim_time=$(grep "Number of generations" $drivers | cut -d':' -f2 | tr -d ' ')
scripts/format_joint_inference_beast_xml.sh ${sim_matrix} ${leaf_tissues} ${template_xml} ${sim_time}

# # run beast joint inference
time java -jar ${metastabayes_jar} -overwrite -working ${dir}/joint_inference_beast.xml > ${dir}/joint_inference_beast_terminal_time.log
beast_posterior_trees="${dir}/joint_inference_beast_tissues.trees"
mcc_tree="${dir}/joint_inference_beast_tissues.tree"
${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${beast_posterior_trees} ${mcc_tree}

# get tissue labeled true tree
conda activate compare_trees
python scripts/format_add_tissues_to_newick.py ${true_tree} ${true_tissues}
true_tissue_tree=${dir}/*_tissue_labeled_tree.nwk

# calculate migration graph F1 scores compared to the true migration graph for machina single result
machina_f1=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${machina_tree} | awk -F' ' '{print $3}')

# calculate migration graph F1 scores compared to the true migration graph for BEAST joint inference MCC single result
python scripts/format_treeannotator_nexus_to_newick.py ${mcc_tree}
beast_mcc_f1=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${mcc_tree}.nwk | awk -F' ' '{print $3}')

# calculate the same F1 score but for sampling all trees from the beast posterior with F1 score weighted by posterior probability
beast_posterior_f1=$(python scripts/migration_graph_f1_true_beast_posterior_trees.py ${true_tissue_tree} ${beast_posterior_trees} | awk -F' ' '{print $3}')

echo "machina"
echo $machina_f1
echo "beast mcc"
echo $beast_mcc_f1
echo "beast posterior"
echo $beast_posterior_f1

echo "${dir},${machina_f1},${beast_mcc_f1}, ${beast_posterior_f1}" >> ${accuracy_file}

conda deactivate
