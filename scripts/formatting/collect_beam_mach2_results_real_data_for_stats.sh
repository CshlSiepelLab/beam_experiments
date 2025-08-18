

dataset_dir_name="serio_prostate_cancer_data"
dataset_name="serio"

indir="/grid/siepel/home/staklins/stored_results/beam/latest_results/$dataset_dir_name/beam_gtr"
outfile="/grid/siepel/home/staklins/stored_results/beam/latest_results/general_graph_stats_for_beam_paper_from_latest_runs_8_2_25/beam_all_results_8_18_25.csv"

# echo "dataset_name,mouse,cp,source_target_edgenum,probability" > "$outfile"

for file in $(find $indir -type f -name "posterior_prob_graph.csv"); do
    mouse=$(basename $(dirname $(dirname $file)))
    cp=$(basename $(dirname $file))
    while IFS= read -r line; do
        echo "$dataset_name,$mouse,$cp,$line" >> "$outfile"
    done < "$file"
done