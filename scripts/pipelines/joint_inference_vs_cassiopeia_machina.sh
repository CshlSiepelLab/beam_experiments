#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

### This pipeline takes in simulated data in the form of an indel character matrix and ground truth tree with tissue labels and the compares cassiopeia->machina and joint tree and tissue BEAST inference method for performance in inferring the migration graph vs the ground truth

# # user inputs
# sim_matrix=$1
# true_tree=$2
# true_tissues=$3

# for testing
sim_matrix="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/987/987_indel_character_matrix.tsv"
true_tree="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/987/tree_seed2064983427.nwk"
true_tissues="sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/987/tree_seed2064983427.vertex.labeling"

# get working dir
dir=$(dirname "$sim_matrix")

# run cassiopeia-greedy on matrix
conda activate simulate
python scripts/cassiopeia_greedy.py $sim_matrix

# run machina on cassiopeia-greedy inferred tree
# Prep cas tree for MACHINA input files
cas_tree_tissues="${dir}/cassiopeia_greedy_inferred.nwk"
machina_dir="${dir}/machina"
primary_tissue="P"
mkdir ${machina_dir}
python ./scripts/machina/prep_machina.py ${cas_tree_tissues} ${machina_dir} ${primary_tissue}
conda deactivate

# Run MACHINA
conda activate machina
#sed -i '0,/0/s/0/GL/' ${machina_dir}/*.tree     # Prevents MACHINA segmentation fault due to input formatting
./scripts/machina/run_machina.sh --edges ${machina_dir}/*.tree --labels ${machina_dir}/*.labeling --colors ${machina_dir}/*_colors.txt --primary-tissue ${primary_tissue} --outdir ${machina_dir}
conda deactivate
# Condense MACHINA output into a labeled tree newick format
conda activate simulate
python ./scripts/machina/post_machina_to_tree.py ${cas_tree_tissues} ${machina_dir}/T-P-0.labeling ${machina_dir}
conda deactivate
# Remove intermediate MACHINA output files
mv ${machina_dir}/machina_tree_all_tissue_labels.nwk ${dir}/
#rm -r ${machina_dir}