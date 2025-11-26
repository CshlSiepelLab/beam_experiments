

indir="/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/migration_windows_sim_data_mig1e-4_11_20_25"

for fitch_tree in $(find $indir/random_consensus_parsimony_tissue_inference -type f -name "parsimony_tissues.nwk"); do
    simname=$(basename $(dirname $fitch_tree))
    echo $simname

    laml_tree=${indir}/laml/${simname}/${simname}_laml_trees.nwk

    python /grid/siepel/home/staklins/projects/crispr_barcode/beam/bayesian_phylogenetic_metastasis/scripts/formatting/transfer_tissue_labels_fitch_to_laml.py \
    $fitch_tree \
    $laml_tree \
    "P" \
    "${indir}/random_consensus_parsimony_tissue_inference/${simname}/parsimony_tissues_with_branch_lengths.nwk"

done


for random_tree in $(find $indir/random_consensus_parsimony_tissue_inference -type f -name "random_tissues.nwk"); do
    simname=$(basename $(dirname $random_tree))
    echo $simname

    laml_tree=${indir}/laml/${simname}/${simname}_laml_trees.nwk

    python /grid/siepel/home/staklins/projects/crispr_barcode/beam/bayesian_phylogenetic_metastasis/scripts/formatting/transfer_tissue_labels_fitch_to_laml.py \
    $random_tree \
    $laml_tree \
    "P" \
    "${indir}/random_consensus_parsimony_tissue_inference/${simname}/random_tissues_with_branch_lengths.nwk"

done


for consensus_tree in $(find $indir/random_consensus_parsimony_tissue_inference -type f -name "consensus_tissues.nwk"); do
    simname=$(basename $(dirname $consensus_tree))
    echo $simname

    laml_tree=${indir}/laml/${simname}/${simname}_laml_trees.nwk

    python /grid/siepel/home/staklins/projects/crispr_barcode/beam/bayesian_phylogenetic_metastasis/scripts/formatting/transfer_tissue_labels_fitch_to_laml.py \
    $consensus_tree \
    $laml_tree \
    "P" \
    "${indir}/random_consensus_parsimony_tissue_inference/${simname}/consensus_tissues_with_branch_lengths.nwk"

done

