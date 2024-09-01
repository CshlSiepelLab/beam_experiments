#!/bin/bash

# required for environment setup on Elzar cshl hpc
source ~/miniconda3/etc/profile.d/conda.sh

### This pipeline takes in simulated data in the form of an indel character matrix and ground truth tree with tissue labels and the compares cassiopeia->machina and joint tree and tissue BEAST inference method for performance in inferring the migration graph vs the ground truth

# user inputs
directory=$1
accuracy_file=$2

# directory="individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24/mS/10965"
# accuracy_file="individual_vs_proper_joint_inference_vs_cassiopeia_machina_6_7_24/individual_joint_inference_accuracy.csv"

sim_matrix=${directory}/*_indel_character_matrix.tsv
true_tree=${directory}/cell_tree_seed*[0-9].nwk
true_tissues=${directory}/cell_tree_seed*.vertex.labeling
leaf_tissues=$(ls ${directory}/cell_tree_*[0-9].labeling)
drivers=${directory}/drivers_seed*.txt
true_migration_graph=${directory}/migration_graph_seed*.csv


# for testing
# sim_matrix="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/24874_indel_character_matrix.tsv"
# true_tree="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874.nwk"
# true_tissues="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874.vertex.labeling"
# leaf_tissues="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/24874/tree_seed24874.labeling"

# get executable paths for beast
treeannotator_path=$(which treeannotator)
logcombiner_path=$(which logcombiner)
metastabayes_jar="../metastabayes/metastabayes.jar"

# get tissue labeled true tree
conda activate compare_trees
python scripts/format_add_tissues_to_newick.py ${true_tree} ${true_tissues}
true_tissue_tree=${directory}/*_tissue_labeled_tree.nwk
conda deactivate

# get working dir
dir=$(dirname "$sim_matrix")

# run cassiopeia-greedy on matrix
conda activate simulate
python scripts/cassiopeia/cassiopeia_greedy.py $sim_matrix

# run machina on cassiopeia-greedy inferred tree
# Prep cas tree for MACHINA input files
cas_tree="${dir}/cassiopeia_greedy_inferred.nwk"
machina_dir="${dir}/machina"
primary_tissue="P"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${cas_tree} ${machina_dir} ${primary_tissue} ${leaf_tissues}
conda deactivate

# Run MACHINA
module load EBModules
module load Gurobi
conda activate machina

#sed -i '0,/0/s/0/GL/' ${machina_dir}/*.tree     # Prevents MACHINA segmentation fault due to input formatting
./scripts/machina/run_machina_tr.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue ${primary_tissue} --outdir ${machina_dir} || echo "Error: Machina execution timed out for ${dir}"

conda deactivate
module unload Gurobi
module unload EBModules

# Condense MACHINA output into a labeled tree newick format
conda activate simulate
python ./scripts/machina/post_machina_tr_to_tree.py ${machina_dir}/P-T-P-R.tree ${machina_dir}/P-T-P-R.labeling ${machina_dir}
conda deactivate
# Remove intermediate MACHINA output files
machina_tree="${dir}/machina_tree_all_tissue_labels.nwk" 
mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${machina_tree}

# setup joint beast inference
conda activate compare_trees
template_xml="inputs/joint_inference_beast_template.xml"
sim_time=$(grep "Number of generations" $drivers | cut -d':' -f2 | tr -d ' ')
beast_dir="${dir}/beast"
mkdir $beast_dir
scripts/format_joint_inference_beast_xml.sh ${sim_matrix} ${leaf_tissues} ${template_xml} ${sim_time} ${beast_dir}

# # run beast joint inference
beast_logs=()
num_chains=5
for ((i=1; i<=$num_chains; i++))
do
  beast_log="${beast_dir}/joint_inference_beast_terminal_time_${i}.log"
  beast_logs+=("$beast_log")
  iter_xml="${beast_dir}/joint_inference_beast_${i}.xml"
  main_xml="${beast_dir}/joint_inference_beast.xml"
  cp $main_xml $iter_xml
  time java -Xmx5g -jar ${metastabayes_jar} -overwrite -working $iter_xml > $beast_log &
done

# Allow for all chains to finish before continuing
wait

# Combine all log and tissue trees files
log_files=""
trees_files=""
ess_convergences=()
for ((i=1; i<=$num_chains; i++))
do
ess_convergence=$(awk '/Operator/ { found=1; next } { if (!found) print }' "${beast_dir}/joint_inference_beast_terminal_time_${i}.log" | awk '{if (NF > 0) print}' | tail -n 1 | awk '{print int($3 + 0.5)}')
ess_convergences+=($ess_convergence)
# only use mcmc chains that have converged since many appear to get lost in the search space
if [[ $ess_convergence -gt 200 ]]; then
  log_files+="-log ${beast_dir}/joint_inference_beast_${i}.log "
  trees_files+="-log ${beast_dir}/joint_inference_beast_${i}_tissues.trees "
fi
done
combined_log="${beast_dir}/joint_inference_beast_combined.log"
combined_trees="${beast_dir}/joint_inference_beast_combined_tissues.trees"
$logcombiner_path $log_files -o $combined_log
$logcombiner_path $trees_files -o $combined_trees

# move combined results to main dir and remove independent chain results
mv $main_xml $dir/
mv $combined_log $dir/
beast_posterior_trees="${dir}/joint_inference_beast_combined_tissues.trees"
mv $combined_trees $beast_posterior_trees

mcc_tree=$(echo "$beast_posterior_trees" | sed 's/.trees/.tree/')
${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${beast_posterior_trees} ${mcc_tree}

# get random and consensus tissue labeled trees
python scripts/consensus_random_tissue_trees.py ${cas_tree} ${leaf_tissues}
random_tissue_tree=${dir}/*_random_tissues.nwk
consensus_tissue_tree=${dir}/*_consensus_tissues.nwk

### F1 scores for the downsampled true cell tree
# calculate migration graph F1 scores compared to the true migration graph for machina single result
tree_machina=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${machina_tree})
tree_machina_f1=$(echo $tree_machina | awk -F' ' '{print $3}')
tree_machina_precision=$(echo $tree_machina | awk -F' ' '{print $5}')
tree_machina_recall=$(echo $tree_machina | awk -F' ' '{print $7}')

# calculate migration graph F1 scores compared to the true migration graph for BEAST joint inference MCC single result
# python scripts/format_treeannotator_nexus_to_newick.py ${mcc_tree}
tree_beast_mcc=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${mcc_tree}.nwk)
tree_beast_mcc_f1=$(echo $tree_beast_mcc | awk -F' ' '{print $3}')
tree_beast_mcc_precision=$(echo $tree_beast_mcc | awk -F' ' '{print $5}')
tree_beast_mcc_recall=$(echo $tree_beast_mcc | awk -F' ' '{print $7}')

# calculate the same F1 score but for sampling all trees from the beast posterior with F1 score weighted by posterior probability
tree_beast_posterior=$(python scripts/migration_graph_f1_true_beast_posterior_trees.py ${true_tissue_tree} ${beast_posterior_trees})
tree_beast_posterior_f1=$(echo $tree_beast_posterior | awk -F' ' '{print $3}')
tree_beast_posterior_precision=$(echo $tree_beast_posterior | awk -F' ' '{print $5}')
tree_beast_posterior_recall=$(echo $tree_beast_posterior | awk -F' ' '{print $7}')
tree_beast_posterior_95ci_binary=$(echo $tree_beast_posterior | awk -F' ' '{print $NF}')

# calculate F1 scores for random and consensus trees
tree_random=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${random_tissue_tree})
tree_random_f1=$(echo $tree_random | awk -F' ' '{print $3}')
tree_random_precision=$(echo $tree_random | awk -F' ' '{print $5}')
tree_random_recall=$(echo $tree_random | awk -F' ' '{print $7}')

tree_consensus=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${consensus_tissue_tree})
tree_consensus_f1=$(echo $tree_consensus | awk -F' ' '{print $3}')
tree_consensus_precision=$(echo $tree_consensus | awk -F' ' '{print $5}')
tree_consensus_recall=$(echo $tree_consensus | awk -F' ' '{print $7}')

### F1 scores for the true migration graph from the entire simulation
true_machina=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_migration_graph} ${machina_tree})
true_machina_f1=$(echo $true_machina | awk -F' ' '{print $3}')
true_machina_precision=$(echo $true_machina | awk -F' ' '{print $5}')
true_machina_recall=$(echo $true_machina | awk -F' ' '{print $7}')

true_beast_mcc=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_migration_graph} ${mcc_tree}.nwk)
true_beast_mcc_f1=$(echo $true_beast_mcc | awk -F' ' '{print $3}')
true_beast_mcc_precision=$(echo $true_beast_mcc | awk -F' ' '{print $5}')
true_beast_mcc_recall=$(echo $true_beast_mcc | awk -F' ' '{print $7}')

true_beast_posterior=$(python scripts/migration_graph_f1_true_beast_posterior_trees.py ${true_migration_graph} ${beast_posterior_trees})
true_beast_posterior_f1=$(echo $true_beast_posterior | awk -F' ' '{print $3}')
true_beast_posterior_precision=$(echo $true_beast_posterior | awk -F' ' '{print $5}')
true_beast_posterior_recall=$(echo $true_beast_posterior | awk -F' ' '{print $7}')
true_beast_posterior_95ci_binary=$(echo $true_beast_posterior | awk -F' ' '{print $NF}')

true_random=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_migration_graph} ${random_tissue_tree})
true_random_f1=$(echo $true_random | awk -F' ' '{print $3}')
true_random_precision=$(echo $true_random | awk -F' ' '{print $5}')
true_random_recall=$(echo $true_random | awk -F' ' '{print $7}')

true_consensus=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_migration_graph} ${consensus_tissue_tree})
true_consensus_f1=$(echo $true_consensus | awk -F' ' '{print $3}')
true_consensus_precision=$(echo $true_consensus | awk -F' ' '{print $5}')
true_consensus_recall=$(echo $true_consensus | awk -F' ' '{print $7}')


# get other stats from the sim
migration_count=$(python scripts/migration_count_from_tree.py $true_tissue_tree | grep -oP 'Migration count: \K.*')
no_nodes_cas_tree=${cas_tree//.nwk/_no_nodes.nwk}
sed 's/node[0-9]*//g' $cas_tree > $no_nodes_cas_tree
cas_rf_dist=$(ete3 compare -t $no_nodes_cas_tree -r $true_tree --unrooted | grep "(..)" | cut -d\| -f 4 | tr -d '[:space:]')
mcc_nwk=${mcc_tree//.tree/.nwk}
joint_rf_dist=$(python scripts/nexus_to_newick.py $mcc_tree | ete3 compare -t $mcc_nwk -r $true_tree --unrooted | grep "(..)" | cut -d\| -f 4 | tr -d '[:space:]')
shannon_mut_matrix=$(python scripts/shannon_entropy_mutation_matrix.py $sim_matrix | grep -oP 'Shannon Entropy scipy: \K.*')
conda deactivate

echo "${dir},\
${tree_machina_f1},\
${tree_machina_precision},\
${tree_machina_recall},\
${tree_beast_mcc_f1},\
${tree_beast_mcc_precision},\
${tree_beast_mcc_recall},\
${tree_beast_posterior_f1},\
${tree_beast_posterior_precision},\
${tree_beast_posterior_recall},\
${tree_random_f1},\
${tree_random_precision},\
${tree_random_recall},\
${tree_consensus_f1},\
${tree_consensus_precision},\
${tree_consensus_recall},\
${tree_beast_posterior_95ci_binary},\
${true_machina_f1},\
${true_machina_precision},\
${true_machina_recall},\
${true_beast_mcc_f1},\
${true_beast_mcc_precision},\
${true_beast_mcc_recall},\
${true_beast_posterior_f1},\
${true_beast_posterior_precision},\
${true_beast_posterior_recall},\
${true_random_f1},\
${true_random_precision},\
${true_random_recall},\
${true_consensus_f1},\
${true_consensus_precision},\
${true_consensus_recall},\
${true_beast_posterior_95ci_binary},\
${ess_convergences[@]},\
${migration_count},\
${cas_rf_dist},\
${joint_rf_dist},\
${shannon_mut_matrix}" >> ${accuracy_file}

# optional clean up of temporary files
# rm -r ${machina_dir}
# rm -r $beast_dir
