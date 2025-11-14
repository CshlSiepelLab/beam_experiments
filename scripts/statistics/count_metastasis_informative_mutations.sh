

indir="/grid/siepel/home/staklins/projects/crispr_barcode/data/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_22_24"
# indir="/grid/siepel/home/staklins/projects/crispr_barcode/data/variable_migration_and_mutation_rates_8_19_24"

script="/grid/siepel/home/staklins/projects/crispr_barcode/beam/bayesian_phylogenetic_metastasis/scripts/statistics/count_metastasis_informative_mutations.py"

for file in $(find $indir -type f -name "cell_tree*.nwk"); do
    vertex_labeling=${file%.nwk}.vertex.labeling
    simname=$(basename $(dirname $file))
    indel_matrix=$(dirname $file)/${simname}_indel_character_matrix.tsv
    outfile=${file%.nwk}_mig_informative_mutation_counts.tsv

    echo ""
    echo "Processing file: $file"
    echo "Vertex labeling: $vertex_labeling"
    echo "Indel matrix: $indel_matrix"
    echo "Output file: $outfile"

    python $script $file $vertex_labeling $indel_matrix $outfile
done


main_outputfile=${indir}/all_mig_informative_mutation_counts.tsv
rm -f $main_outputfile
files=$(find $indir -type f -name "*_mig_informative_mutation_counts.tsv")
echo -e "simname\t$(head -n 1 $(echo $files | awk '{print $1}'))" > $main_outputfile
for f in $files; do
    simname=$(basename $(dirname $f))
    tail -n +2 $f | awk -v simname=$simname '{print simname"\t"$0}' >> $main_outputfile
done
