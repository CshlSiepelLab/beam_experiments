#!/bin/bash

### This script is to simulate trees with node tissue labels according to the modified TTP / MACHINA agent based model of tumor growth and metastasis

# set path to c++ executable for agent based model of cancer growth and metastasis
simulate_executable="/grid/siepel/home_norepl/staklins/barcode_simulator/scripts/simulator/build/simulate"

# set dir to hold all sims
sim_dir="sim_trees_3_25_24"
mkdir ${sim_dir}

# make pattern directories
mkdir ${sim_dir}/mS
mkdir ${sim_dir}/pS
mkdir ${sim_dir}/pM
mkdir ${sim_dir}/pR

# set migration pattern options
migration_patterns=(0 1 2 3)
for pattern in ${migration_patterns[@]}; do
for ((i = 0; i < 10; i++)); do

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
num_cells=1000
max_anatomical_sites=10
migration_rate="1e-3"   ### amped up migration rate from the default for testing purposes of model fit
carrying_capacity="5e2" ### reduced carying capacity to a similar degree as migration rate to model smaller cancer population with similar dynamics as a larger one

# save running parameters to file
run_conditions_file="${outprefix}/sim_run_conditions.txt"
echo -e "sim_name\t${sim_dir}\n
output_directory\t${outprefix}\n
seed\t${seed}\n
pattern\t${pattern}\n
num_cells\t${num_cells}\n
migration_rate\t${migration_rate}\n
carrying_capacity\t${carrying_capacity}\n
max_anatomical_sites\t${max_anatomical_sites}" > ${run_conditions_file}

# run simulator
${simulate_executable} -C ${num_cells} -K ${carrying_capacity} -m ${max_anatomical_sites} -p ${pattern} -mig ${migration_rate} -s ${seed} -o ${outprefix}

# prep machina output
tree_file="${outprefix}/*.tree"
label_file="${outprefix}/*.vertex.labeling"
sed -i 's/\//;/g' ${tree_file}
sed -i 's/\//;/g' ${label_file}
python ./scripts/machina_sims_to_newick_format.py ${tree_file} ${label_file}

done
done