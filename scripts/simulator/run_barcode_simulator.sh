
barcode_simulator_dir="../barcode_simulator/scripts/simulator"

# set dir to hold all sims
sim_dir="data/test_150_tips_9_12_24"
mkdir -p ${sim_dir}

migration_rates=(1e-6)
# migration_rates=(1e-4 1e-5 1e-6 1e-7)

mutation_rates=(0.005)
# mutation_rates=(0.001 0.005 0.01)

for migration_rate in ${migration_rates[@]}; do

# migration_rate="1e-6,0,1e-6,1e-6,1e-6"

# m=$(echo $migration_rate | sed 's/1e-//g')

for mutrate in ${mutation_rates[@]}; do

# j=$(echo $mutrate | sed 's/^0\.//')

for ((i = 0; i < 5; i++)); do

seed=$RANDOM


outprefix="${HOME}/bayesian_phylogenetic_metastasis/${sim_dir}/${seed}"
# outprefix="${HOME}/bayesian_phylogenetic_metastasis/${sim_dir}/mig${m}_mut${j}_${seed}"

# re-draw seed if it already exists
while [ -d "${outprefix}" ]; do
    seed=$RANDOM
    outprefix="${HOME}/bayesian_phylogenetic_metastasis/${sim_dir}/${seed}"
    # outprefix="${HOME}/bayesian_phylogenetic_metastasis/${sim_dir}/mig${m}_mut${j}_${seed}"
done

mkdir ${outprefix}

# set simulator parameters
num_generations=250
max_anatomical_sites=-1
# migration_rate="1e-6"
mutFreqThreshold=0.05
carryingCapacity="5e4"
driverProb="1e-7"
num_cells_downsample=100

num_sites=50
design="RANDOM"

# save running parameters to file
run_conditions_file="${outprefix}/sim_run_conditions.txt"
touch $run_conditions_file
echo -e "sim_name\t${sim_dir}\n
output_directory\t${outprefix}\n
seed\t${seed}\n
pattern\t${pattern}\n
num_cells\t${num_cells}\n
migration_rate\t${migration_rate}\n
num_sites\t${num_sites}\n
mutrate\t${mutrate}\n
num_cells_downsample\t${num_cells_downsample}\n
transition_probs_file\t${transition_probs}" > ${run_conditions_file}

# run simulator
echo "Running barcode simulator with seed ${seed} and migration pattern ${pattern}"
echo "$barcode_simulator_dir/build/simulate -C ${num_generations} -mig ${migration_rate} -s ${seed} -o ${outprefix} -m ${max_anatomical_sites} -f ${mutFreqThreshold} -K ${carryingCapacity} -D ${driverProb} -d ${num_cells_downsample}" >> ${run_conditions_file}
timeout 10m $barcode_simulator_dir/build/simulate -C ${num_generations} -mig ${migration_rate} -s ${seed} -o ${outprefix} -m ${max_anatomical_sites} -f ${mutFreqThreshold} -K ${carryingCapacity} -D ${driverProb} -d ${num_cells_downsample} >> ${run_conditions_file} || (rm -rf ${outprefix} && continue)

# if no migrations, then repeat simulation
migrations=$(grep -v '^$' ${outprefix}/migration_graph* | wc -l)
if [ $migrations -lt 2 ]; then
    echo "No migrations occurred, repeating simulation"
    rm -rf ${outprefix}
    i=$((i-1))
    continue
fi


# run barcode simulator to overlay barcode data for machina simulator tree and tissues output
machina_tree=${outprefix}/cell_tree_seed*.nwk
machina_tissues=${outprefix}/cell_tree_seed*.vertex.labeling
echo "Running overlay barcode data for machina simulator tree and tissues output"
echo "python $barcode_simulator_dir/overlay_barcode_machina_simulator.py ${outprefix} ${design} ${num_sites} ${mutrate} ${machina_tree} ${machina_tissues}" >> ${run_conditions_file}
python $barcode_simulator_dir/overlay_barcode_machina_simulator.py ${outprefix} ${design} ${num_sites} ${mutrate} ${machina_tree} ${machina_tissues}

done
done
done
