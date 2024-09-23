#!/bin/bash

# use matplotlib conda env activated when running this script

true_trees=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_9_12_24/raw_data -type f -name "*issue_labeled_tree.nwk")
true_trees2=$(echo $true_trees | sed 's/ /,/g')

python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/plotting/plot_precision_recall_downsampling.py \
$true_trees2 \
P \
/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_150_tips_9_12_24/precision_recall \
