#!/bin/bash

# need an environment with ete3 installed
# mamba activate compare_trees

posterior_files=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/snakemake_performance_repeat_origin_scaling_implemented_10_15_24_uniform_50cells_50sites_data_7_24_24/metastabayes -type f -name "combined.trees")

primary_tissue="P"

for posterior_file in $posterior_files; do
    echo "Starting $posterior_file"
    dir=$(dirname $posterior_file)
    working_dir=${dir}/posterior_expected_over_parsimony
    mkdir -p $working_dir
    echo -e "id,beast_migration_count,parsimony_migration_count,posterior_excess" > $working_dir/posterior_expected_migration_counts_over_parsimony.csv

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
        beast_migration_count=$(python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/migration_count_from_tree.py $beast_result $primary_tissue | cut -d' ' -f3)
        parsimony_migration_count=$(python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/statistics/migration_count_from_tree.py $parsimony_result $primary_tissue | cut -d' ' -f3)
        excess=$((beast_migration_count - parsimony_migration_count))

        # save results
        echo -e "$id,$beast_migration_count,$parsimony_migration_count,$excess" >> $working_dir/posterior_expected_migration_counts_over_parsimony.csv

        # clean up
        rm -r $working_dir_id
    done < <(grep 'tree STATE' "$posterior_file")
done