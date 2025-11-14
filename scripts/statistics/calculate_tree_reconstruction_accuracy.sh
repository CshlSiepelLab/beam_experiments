
# # For ideal sims
# data_indir="/grid/siepel/home/staklins/projects/crispr_barcode/data/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_22_24"
# results_indir="/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24"

# For variable sims
data_indir="/grid/siepel/home/staklins/projects/crispr_barcode/data/variable_migration_and_mutation_rates_8_19_24"
results_indir="/grid/siepel/home/staklins/projects/crispr_barcode/results/beam/latest_results/variable_migration_and_mutation_rates_data_8_19_24"


script="/grid/siepel/home/staklins/projects/crispr_barcode/beam/bayesian_phylogenetic_metastasis/scripts/statistics/calculate_tree_reconstruction_accuracy.py"

outdir=${results_indir}/tree_reconstruction_accuracy
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

main_outfile=${outdir}/all_tree_reconstruction_accuracy.csv
if [ -f $main_outfile ]; then
    rm $main_outfile
fi

for true_tree_file in $(find $data_indir -type f -name "cell_tree*.nwk"); do

    simname=$(basename $(dirname $true_tree_file))

    laml_tree_file=${results_indir}/laml/${simname}/${simname}_laml_trees.nwk
    cassiopeia_tree_file=${results_indir}/cassiopeia_greedy/${simname}/cassiopeia_greedy_inferred.nwk
    beam_trees_file=${results_indir}/beam_gtr/${simname}/combined.trees

    outfile=${outdir}/${simname}_tree_reconstruction_accuracy.csv

    pathfinder_tree_file=${results_indir}/pathfinder/${simname}/scratc*/clone_aln.nwk
    if [ ! -f $pathfinder_tree_file ]; then
        pathfinder_tree_file=None
    fi

    python $script $true_tree_file $laml_tree_file $cassiopeia_tree_file $beam_trees_file $outfile $pathfinder_tree_file

    if [ ! -f $main_outfile ]; then
        echo "simname,$(head -n 1 $outfile)" > $main_outfile
    fi
    tail -n +2 $outfile | awk -v simname="$simname" '{print simname "," $0}' >> $main_outfile

    rm $outfile
done


