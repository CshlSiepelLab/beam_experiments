
barcode_simulator_dir="/grid/siepel/home/staklins/projects/crispr_barcode/barcode_simulator/scripts/simulator"

# set dir to hold all sims
outdir="/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_22_24"
mkdir -p ${outdir}

migration_rates=(1e-6)
mutation_rates=(0.0025)

for migration_rate in ${migration_rates[@]}; do

    for mutrate in ${mutation_rates[@]}; do

        for ((i = 0; i < 80; i++)); do

            echo "Running simulation $i"

            # Get a random seed that does not already exist
            seed=$RANDOM
            outprefix="${outdir}/${seed}"
            while [ -d "${outprefix}" ]; do
                seed=$RANDOM
                outprefix="${outdir}/${seed}"
            done
            mkdir ${outprefix}
            log_file="${outprefix}/log.txt"

            # set simulator parameters
            num_generations=250
            max_anatomical_sites=-1
            num_cells_downsample=50
            num_sites=50

            # run simulator
            timeout 10m $barcode_simulator_dir/build/simulate -c ${num_generations} -mig ${migration_rate} -s ${seed} -o ${outprefix} -m ${max_anatomical_sites} -d ${num_cells_downsample} >> ${log_file} || (rm -rf ${outprefix} && continue)

            # if no migrations, then repeat simulation
            migrations=$(grep -v '^$' ${outprefix}/migration_graph* | wc -l)
            if [ $migrations -lt 2 ]; then
                echo "No migrations occurred, repeating simulation"
                rm -rf ${outprefix}
                i=$((i-1))
                continue
            fi


            # run barcode simulator to overlay barcode data
            python $barcode_simulator_dir/overlay_barcode.py ${outprefix} ${num_sites} ${mutrate} ${outprefix}/cell_tree_seed*.nwk

        done
    done
done
