#!/bin/bash

### Basic pipeline to take the posterior files of Beast joint tree and migration inference and run treeannotator to calculate true tree vs MCC tree F1 score for the migration graph

input_dir="compare_mcmc_length_joint_inference_4_26_24"
input_posterior_files=$(find $input_dir -type f -name "*tissues.trees")
# input_posterior_files=(compare_mcmc_length_joint_inference_4_26_24/joint_inference_beast_100000_tissues.trees compare_mcmc_length_joint_inference_4_26_24/joint_inference_beast_1million_tissues.trees compare_mcmc_length_joint_inference_4_26_24/joint_inference_beast_10million_tissues.trees compare_mcmc_length_joint_inference_4_26_24/joint_inference_beast_25million_tissues.trees)

# find treeannotator to condense beast posterior to mcc
treeannotator=$(which treeannotator)

# specify true tree with tissue labels for comparisons, with one tree for all runs in this case but could make this seperate true trees to match runs by doing loop with an iterator and indexing 2 arrays
true_tree="/home/staklins/bayesian_phylogenetic_metastasis/joint_inference_vs_cassiopeia_machina_cellTree_simdataset_4_25_24/pS/4011/cell_tree_seed4011_tissue_labeled_tree.nwk"

for file in ${input_posterior_files[@]}
do
# get mcc tree
mcc_tree=$(echo $file | sed 's/.trees/_mcc.tree/g')
$treeannotator -burnin 10 -topology MCC -height mean -file $file $mcc_tree

# format mcc to newick
python scripts/format_treeannotator_nexus_to_newick.py ${mcc_tree}

# get mcc f1 vs true
beast_mcc_f1=$(python scripts/migration_graph_f1_true_inferred_trees.py $true_tree $mcc_tree.nwk | awk -F' ' '{print $3}')

# get fractional posterior f1 vs true
beast_posterior_f1=$(python scripts/migration_graph_f1_true_beast_posterior_trees.py $true_tree $file | awk -F' ' '{print $3}')

# output f1 scores to file
outputfile=$(echo $file | sed 's/.trees/_f1scores.txt/g')
echo "$outputfile"
echo -e "file,mccf1,posteriorf1\n$file,$beast_mcc_f1,$beast_posterior_f1" > $outputfile
done