
main_dir="/grid/siepel/home/staklins/stored_results/beam/latest_results/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_24_24"


primary_tissue="P"
outfile="$main_dir/true_tree_stats.txt"

echo -e "sim_name,migration_count,comigration_count,num_multiedges,met_to_met,reseeding,clonality" > $outfile

# Get all of the true tree labeled newick files to classify the true migration graphs
files=$(find $main_dir/true_trees -type f -name "tissue_labeled_tree.nwk")

for file in $files; do

    sim_name=$(dirname $file | rev | cut -d'/' -f1 | rev)

    # get stats
    python ./scripts/formatting/call_migration_comigration_multiedges_topology_from_labeled_newick.py $sim_name $file $primary_tissue $outfile
done

