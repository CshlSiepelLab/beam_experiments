#!/bin/bash

### This script is to simulate trees with node tissue labels according to the modified TTP / MACHINA agent based model of tumor growth and metastasis

# set path to c++ executable for agent based model of cancer growth and metastasis
simulate_executable="/grid/siepel/home_norepl/staklins/machina/build/simulate"

# set dir to hold all sims
sim_dir="sim_trees_originalMACHINA_3_27_24"
mkdir ${sim_dir}

# make pattern directories
mkdir ${sim_dir}/mS
mkdir ${sim_dir}/pS
mkdir ${sim_dir}/pM
mkdir ${sim_dir}/pR

# set migration pattern options
migration_patterns=(0 1 2 3)
for pattern in ${migration_patterns[@]}; do
for ((i = 0; i < 20; i++)); do

# make dir specific to the seed number for each sim and the migration pattern
seed=$RANDOM
if (( $pattern == 0 )); then
    pattern_dir="mS"
elif (( $pattern == 1 )); then
    pattern_dir="pS"
elif (( $pattern == 2 )); then
    pattern_dir="pM"
elif (( $pattern == 3 )); then
    pattern_dir="pR"
fi

outprefix="${HOME}/bayesian_phylogenetic_metastasis/${sim_dir}/${pattern_dir}/${seed}"
mkdir ${outprefix}

# set simulator parameters
max_anatomical_sites=8
migration_rate="1e-3"
coloring_file="/grid/siepel/home_norepl/staklins/machina/data/sims/coloring.txt"

# save running parameters to file
run_conditions_file="${outprefix}/sim_run_conditions.txt"
echo -e "sim_name\t${sim_dir}\n
output_directory\t${outprefix}\n
seed\t${seed}\n
pattern\t${pattern}\n
migration_rate\t${migration_rate}\n
max_anatomical_sites\t${max_anatomical_sites}" > ${run_conditions_file}

# run simulator
${simulate_executable} -C 200 -c ${coloring_file} -m ${max_anatomical_sites} -p ${pattern} -mig ${migration_rate} -s ${seed} -o ${outprefix}

# prep machina output
tree_file="${outprefix}/T_*.tree"
label_file="${outprefix}/T_*.vertex.labeling"

python ./scripts/machina_sims_to_newick_format.py ${tree_file} ${label_file}

done
done