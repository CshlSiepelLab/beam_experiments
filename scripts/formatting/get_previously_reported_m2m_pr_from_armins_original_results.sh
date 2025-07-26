

maindir="/grid/siepel/home/staklins/stored_results/cancer_discovery_serio_paper/machina"

processing_script="/grid/siepel/home/staklins/projects/crispr_barcode/bayesian_phylogenetic_metastasis/scripts/formatting/get_previously_reported_m2m_pr_from_armins_original_results.py"

primary_tissue="PRL"

outfile="/grid/siepel/home/staklins/stored_results/bayesian_migration_graph_inference/latest_results/previously_reported_m2m_pr_armins_original_serio_data_7_26_25/all_cp_classifications.csv"
mkdir -p "$(dirname "$outfile")"

echo "mmus,cp,m2m,pr" > "$outfile"

edge_files=$(find "$maindir" -name "G*binarized_parent_child.txt" -type f)

for file in $edge_files; do
    echo $file

    labeling_file="${file/parent_child.txt/labels.txt}"
    echo $labeling_file

    result=$(python $processing_script $file $labeling_file $primary_tissue)

    m2m=$(echo $result | cut -d' ' -f4)
    pr=$(echo $result | cut -d' ' -f2)

    mmus=$(basename $(dirname $(dirname $(dirname $file))) | cut -d'_' -f1)
    cp=$(basename $(dirname $file) | cut -d'_' -f1)

    echo "$mmus,$cp,$m2m,$pr" >> "$outfile"
done