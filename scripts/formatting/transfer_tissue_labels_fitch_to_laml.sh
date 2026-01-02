

indir="/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/migration_windows_sim_data_mig1e-6_150tips_12_30_25"

script="/grid/siepel/home/staklins/projects/crispr_barcode/beam/bayesian_phylogenetic_metastasis/scripts/formatting/transfer_tissue_labels_fitch_to_laml.py"

primary_tissue="P"

for fitch_tree in $(find $indir/random_consensus_parsimony_tissue_inference -type f -name "parsimony_tissues.nwk"); do
    simname=$(basename $(dirname $fitch_tree))
    echo $simname

    laml_tree=${indir}/laml/${simname}/${simname}_laml_trees.nwk

    python $script $fitch_tree $laml_tree $primary_tissue "${indir}/random_consensus_parsimony_tissue_inference/${simname}/parsimony_tissues_with_branch_lengths.nwk"
done


for random_tree in $(find $indir/random_consensus_parsimony_tissue_inference -type f -name "random_tissues.nwk"); do
    simname=$(basename $(dirname $random_tree))
    echo $simname

    laml_tree=${indir}/laml/${simname}/${simname}_laml_trees.nwk

    python $script $random_tree $laml_tree $primary_tissue "${indir}/random_consensus_parsimony_tissue_inference/${simname}/random_tissues_with_branch_lengths.nwk"

done


for consensus_tree in $(find $indir/random_consensus_parsimony_tissue_inference -type f -name "consensus_tissues.nwk"); do
    simname=$(basename $(dirname $consensus_tree))
    echo $simname

    laml_tree=${indir}/laml/${simname}/${simname}_laml_trees.nwk

    python $script $consensus_tree $laml_tree $primary_tissue "${indir}/random_consensus_parsimony_tissue_inference/${simname}/consensus_tissues_with_branch_lengths.nwk"
done

