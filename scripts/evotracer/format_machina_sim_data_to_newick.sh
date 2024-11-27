#!/bin/bash
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate ete3

# failed_seeds_file="failed_seeds.txt"
# good_seeds_file="good_seeds.txt"
# echo "" > $failed_seeds_file
# echo "" > $good_seeds_file

for dir in */; do
    dir_prefix=$(echo $dir | awk -F'/' '{print $1}')
    tree_file="${dir}T_${dir_prefix}.tree"
    label_file="${dir}T_${dir_prefix}.vertex.labeling"
    echo $dir
    python ../scripts/machina_sims_to_newick_format.py ${tree_file} ${label_file}
    # Check if the Python command failed (non-zero exit code)
    # if [ $? -ne 0 ]; then
    #     echo $dir >> $failed_seeds_file
    # else
    #     echo $dir >> $good_seeds_file
    # fi
done