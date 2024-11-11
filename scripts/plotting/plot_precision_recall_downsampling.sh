#!/bin/bash

# use matplotlib conda env activated when running this script

true_trees=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_100_tips_10_31_24/raw_data -type f -name "*issue_labeled_tree.nwk")
true_trees2=$(echo $true_trees | sed 's/ /,/g')

# to instead subset to test a single tree
# true_trees2=$(echo "$true_trees" | tr ' ' '\n' | head -n 1 | tr '\n' ',' | sed 's/,$//')

python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/plotting/plot_precision_recall_downsampling.py \
$true_trees2 \
P \
/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/downsampling_test_100_tips_10_31_24/precision_recall \
