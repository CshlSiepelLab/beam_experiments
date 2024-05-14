#!/bin/bash

outfile="results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/rf_distances_true_cassiopeia.csv"

true_tree_files=$(find results/moreSims_joint_inference_vs_cassiopeia_machina_vs_random_cellTree_simdataset_5_3_24/ -type f -name cell_tree_seed*[0-9].nwk)

for file in $true_tree_files; do
# get cassiopeia file
cas_tree_file=$(echo $file | cut -d\/ -f 1-4)/cassiopeia_greedy_inferred.nwk
new_cas_tree_file=$(echo $cas_tree_file | sed 's/\.nwk/_no_nodes.nwk/')
cat $cas_tree_file | sed 's/node[0-9]*//g' > $new_cas_tree_file
echo $file >> $outfile
echo $new_cas_tree_file >> $outfile
ete3 compare -t $new_cas_tree_file -r $file --unrooted >> $outfile
done

