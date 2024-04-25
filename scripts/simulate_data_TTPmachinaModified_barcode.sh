#!/bin/bash

### This script is to simulate trees with node tissue labels according to the modified TTP / MACHINA agent based model of tumor growth and metastasis
### need to activate conda environment names simulate that has cassiopeia and ete3

barcode_simulator_dir="../barcode_simulator/scripts/simulator"

# set dir to hold all sims
sim_dir="sim_data_50cellTrees_barcodes_4_25_24"
mkdir ${sim_dir}

# make pattern directories
mkdir ${sim_dir}/mS
mkdir ${sim_dir}/pS
mkdir ${sim_dir}/pM
mkdir ${sim_dir}/pR

# set migration pattern options
migration_patterns=(0 1 2 3)
#migration_patterns=(0)

mutation_rates=(0.01)

for mutrate in ${mutation_rates[@]}; do
for pattern in ${migration_patterns[@]}; do
for ((i = 0; i < 5; i++)); do

# make dir specific to the seed number for each sim and the migration pattern
if (( $pattern == 0 )); then
    pattern_dir="mS"
elif (( $pattern == 1 )); then
    pattern_dir="pS"
elif (( $pattern == 2 )); then
    pattern_dir="pM"
elif (( $pattern == 3 )); then
    pattern_dir="pR"
fi

seed=$RANDOM

outprefix="${HOME}/bayesian_phylogenetic_metastasis/${sim_dir}/${pattern_dir}/${seed}"
mkdir ${outprefix}

# set simulator parameters
num_cells=-1
max_anatomical_sites=8
migration_rate="1e-6"
mutFreqThreshold=0.05
carryingCapacity="5e4"
driverProb="1e-7"
num_cells_downsample=50

num_sites=10
design="RANDOM"


# save running parameters to file
run_conditions_file="${outprefix}/sim_run_conditions.txt"
echo -e "sim_name\t${sim_dir}\n
output_directory\t${outprefix}\n
seed\t${seed}\n
pattern\t${pattern}\n
num_cells\t${num_cells}\n
migration_rate\t${migration_rate}\n
num_sites\t${num_sites}\n
mutrate\t${mutrate}\n
num_cells_downsample\t${num_cells_downsample}" > ${run_conditions_file}

# run simulator
echo "Running barcode simulator with seed ${seed} and migration pattern ${pattern}"
echo "$barcode_simulator_dir/build/simulate -C ${num_cells} -p ${pattern} -mig ${migration_rate} -s ${seed} -o ${outprefix} -m ${max_anatomical_sites} -f ${mutFreqThreshold} -K ${carryingCapacity} -D ${driverProb} -d ${num_cells_downsample}" >> ${run_conditions_file}
timeout 1m $barcode_simulator_dir/build/simulate -C ${num_cells} -p ${pattern} -mig ${migration_rate} -s ${seed} -o ${outprefix} -m ${max_anatomical_sites} -f ${mutFreqThreshold} -K ${carryingCapacity} -D ${driverProb} -d ${num_cells_downsample} >> ${run_conditions_file} || (rm -rf ${outprefix} && continue)


# run barcode simulator to overlay barcode data for machina simulator tree and tissues output
machina_tree=${outprefix}/cell_tree_seed*.nwk
machina_tissues=${outprefix}/cell_tree_seed*.vertex.labeling
echo "Running overlay barcode data for machina simulator tree and tissues output"
echo "python $barcode_simulator_dir/overlay_barcode_machina_simulator.py ${outprefix} ${design} ${num_sites} ${mutrate} ${machina_tree} ${machina_tissues}" >> ${run_conditions_file}
python $barcode_simulator_dir/overlay_barcode_machina_simulator.py ${outprefix} ${design} ${num_sites} ${mutrate} ${machina_tree} ${machina_tissues}

done
done
done
