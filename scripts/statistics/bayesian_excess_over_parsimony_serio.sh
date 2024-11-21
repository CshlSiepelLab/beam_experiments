#!/bin/bash

# need an environment with ete3 installed
# mamba activate compare_trees

posterior_file="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24/beam/MMUS1457/CP01/combined.trees"
primary_tissue="PRL"
echo "Starting $posterior_file"
dir=$(dirname $posterior_file)
name=$(basename $dir)
working_dir=${dir}/posterior_expected_over_parsimony
mkdir -p $working_dir
echo -e "name,sample_id,beast_migration_count,parsimony_migration_count" > $working_dir/posterior_expected_migration_counts_over_parsimony.csv

while IFS= read -r line; do
    id=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
    working_dir_id=${working_dir}/${id}
    mkdir -p $working_dir_id
    newick=$(echo $line | cut -d' ' -f5-)

    # process beast tree and get migration count
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/formatting/beast_posterior_tree_to_newicks.py $newick $id $working_dir_id
    beast_result=$working_dir_id/${id}_beast.newick
    plain_newick=$working_dir_id/${id}.newick
    tip_tissues=$working_dir_id/${id}_tip_tissues.tsv

    # get parsimony result
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/consensus_random/parsimony_only.py $plain_newick $tip_tissues $working_dir_id $primary_tissue
    parsimony_result=$working_dir_id/parsimony_tissues.nwk

    # get migration counts
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/migration_count_from_tree.py $beast_result $primary_tissue > $working_dir_id/beast_count.txt
    beast_migration_count=$(cat $working_dir_id/beast_count.txt)
    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/migration_count_from_tree.py $parsimony_result $primary_tissue > $working_dir_id/parsimony_count.txt
    parsimony_migration_count=$(cat $working_dir_id/parsimony_count.txt)

    # save results
    echo -e "$name,$id,$beast_migration_count,$parsimony_migration_count" >> $working_dir/posterior_expected_migration_counts_over_parsimony.csv

    # clean up
    rm -r $working_dir_id
done < <(grep 'tree STATE' "$posterior_file")
