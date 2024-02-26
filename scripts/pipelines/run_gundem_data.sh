#!/bin/bash

# Necessary line to access conda commands for bash script on CSHL HPC cluster
# source ~/anaconda3/etc/profile.d/conda.sh

# Necessary line to access conda commands on Evolgen lab server (need to make these the same long term)
source ~/miniconda3/etc/profile.d/conda.sh

pipeline_run="beast_gundem_2015_2_21_24"
mkdir ${pipeline_run}

symmetric=(true false)

for sym in "${symmetric[@]}"; do

for tree_file in machina_data/realdata/gundem_2015/reported_clonetrees/*.tree; do

dir_prefix=$(basename "$tree_file" .tree)
labeling_file="${tree_file%.tree}.labeling"
primary_tissue="prostate"

# Rename dir for running both symmetrical and asymmetrical beast in parallel
if [ "$sym" = true ]; then
    sym_name="sym"
else
    sym_name="asym"
fi
sym_prefix="${dir_prefix}_${sym_name}"

# Make directory specific for each patient
dir="${pipeline_run}/${sym_prefix}"
mkdir ${dir}

# Copy original files into working directory for each patient
cp ${tree_file} ${dir}/
cp ${labeling_file} ${dir}/

tree_file="${dir}/${dir_prefix}.tree"
labeling_file="${tree_file%.tree}.labeling"

conda activate ete3

python scripts/machina_realdata_to_newick.py $tree_file $labeling_file $primary_tissue

tree="${dir}/${dir_prefix}_unlabeled_tree.nwk"
tissues="${dir}/${dir_prefix}_tissues.tsv"

python ./scripts/format_xml_template_inputs_fixedTreeAnalysis_from_sim.py ${tree} ${tissues}

conda deactivate

# Format template xml
seqfile="${dir}/${dir_prefix}_unlabeled_tree_sequences_formatted_for_xml.txt"
taxafile="${dir}/${dir_prefix}_unlabeled_tree_taxonset_formatted_for_xml.txt"
traitfile="${dir}/${dir_prefix}_unlabeled_tree_traitset_formatted_for_xml.txt"
newickfile="${dir}/${dir_prefix}_unlabeled_tree_newick_formatted_for_xml.txt"

xml_template="inputs/template_xml_fixedtreeanalysis_machina_sim_universal.xml"
symmetric="${sym}"

# Shorter chain length for symmetrical setup since convergence is reached earlier with less parameters
if [ "$sym" = true ]; then
    chainlength=1000000
else
    chainlength=10000000
fi

scripts/format_template_fixedTreeAnalysis_xml_from_sim.sh ${seqfile} ${taxafile} ${traitfile} ${newickfile} ${xml_template} ${primary_tissue} ${symmetric} ${chainlength}

# Run BEAST2 on formatted xml with output automatically in sim directory
beast_path=$(which beast)
# start_time=$(date +%s.%N)
${beast_path} -overwrite -working ${dir}/${dir_prefix}_unlabeled_tree_final_input_xml.xml
# end_time=$(date +%s.%N)
# beast_time=$(printf "%.2f" $(echo "$end_time - $start_time" | bc))

# Get Maximum Clade Credibility tree from posterior of trees
treeannotator_path=$(which treeannotator)
consensus_tree="${dir}/tissue_tree_with_trait.tree"
${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${dir}/tissue_tree_with_trait.trees ${consensus_tree}

# Plot rate matrix and FigTree from BEAST results
conda activate ggplot2
logfile="${dir}/${dir_prefix}_unlabeled_tree_final_input_xml.log"
Rscript scripts/plot_rate_matrix_from_beast_log.R $logfile $primary_tissue
conda deactivate

treefile="${dir}/tissue_tree_with_trait.tree"
scripts/plot_figtree_consensus_tree.sh $treefile

conda activate ggtree
Rscript scripts/plot_tree_piecharts_ggtree.R ${consensus_tree} ${primary_tissue}
conda deactivate


### NEED TO SETUP MACHINA BELOW AND CHECK FOR BUGS IF USING THIS TO COMPARE. USING ORIGINAL MACHINA FIGURES FOR NOW TO COMPARE WITH BEAST.
# # Prep and run MACHINA to compare
# # Prep simulated true tree for MACHINA input files
# labeled_tree="${dir}/${dir_prefix}_tissue_labeled_tree.nwk"
# machina_dir="${dir}/machina"
# mkdir ${machina_dir}
# conda activate simulate
# python ./scripts/machina/prep_machina.py ${labeled_tree} ${machina_dir} ${primary_tissue}
# conda deactivate

# # Run MACHINA
# conda activate machina
# sed -i '0,/0/s/0/GL/' ${machina_dir}/*.tree     # Prevents MACHINA segmentation fault due to input formatting
# start_time=$(date +%s.%N)
# ./scripts/machina/run_machina.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue ${primary_tissue} --outdir ${machina_dir}
# end_time=$(date +%s.%N)
# machina_time=$(printf "%.2f" $(echo "$end_time - $start_time" | bc))
# conda deactivate

# # Condense MACHINA output into a labeled tree newick format
# conda activate simulate
# python ./scripts/machina/post_machina_to_tree.py ${labeled_tree} ${machina_dir}/T-P-0.labeling ${machina_dir}
# conda deactivate

# # Remove intermediate MACHINA output files
# mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${dir}/
# rm -r ${machina_dir}

done
done