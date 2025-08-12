
working_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50"
primaryTissue="PRL"

# find parsimony files since they always complete for all methods
parsimony_files=$(find $working_dir/parsimony_tissue_inference -type f -name parsimony_tissues_random.nwk)

process_file() {
    file=$1
    echo $file
    parsimony_dir=$(dirname $file)
    parsimony=$parsimony_dir/all_parsimony_solutions.txt
    if [ -f "$parsimony" ]; then
        rm "$parsimony"
    fi

    all_parsimony_solutions=$(find $parsimony_dir -type f -name "parsimony_tissues_all_solutions*.nwk")
    for solution in $all_parsimony_solutions; do
        cat $solution >> $parsimony
    done

    mouse=$(echo $parsimony_dir | rev | cut -d'/' -f 2 | rev)
    cp=$(echo $parsimony_dir | rev | cut -d'/' -f 1 | rev)
    name=${mouse}_${cp}

    machina=$working_dir/machina/$mouse/$cp/${primaryTissue}-G-${primaryTissue}-R.tree
    metient=$working_dir/metient/$mouse/$cp/${mouse}_${cp}_${primaryTissue}_migration_graphs.txt
    beam=$working_dir/beam/$mouse/$cp/combined.trees
    outdir=$working_dir/migration_count_all_method_comparison
    mkdir -p $outdir

    python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/formatting/call_migration_counts_machina_metient_beam.py \
    $machina \
    $metient \
    $beam \
    $parsimony \
    $primaryTissue \
    $outdir  \
    $name
}

export -f process_file
export working_dir
export primaryTissue

cores=80
echo "$parsimony_files" | parallel -j $cores process_file