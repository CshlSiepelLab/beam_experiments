
barcode_simulator_dir="/grid/siepel/home/staklins/projects/crispr_barcode/barcode_simulator/scripts/simulator"

# set dir to hold all sims
outdir="/grid/siepel/home/staklins/stored_data/crispr_barcode_related_data/uniform_50cells_50sites_0.0025mut_10-6mig_data_8_22_24"
mkdir -p ${outdir}

migration_rates=(1e-6)
mutation_rates=(0.0025)

for migration_rate in ${migration_rates[@]}; do

    for mutrate in ${mutation_rates[@]}; do

        for ((i = 0; i < 1; i++)); do

            echo "Running simulation $i"

            # Get a random seed that does not already exist
            seed=$RANDOM
            outprefix="${outdir}/${seed}"
            while [ -d "${outprefix}" ]; do
                seed=$RANDOM
                outprefix="${outdir}/${seed}"
            done
            mkdir ${outprefix}

            # set simulator parameters
            num_generations=250
            max_anatomical_sites=-1
            num_cells_downsample=50
            num_sites=50

            # run executable from beam_sup package
            run_met_cancer_barcode_simulations.sh -\
            -outdir ${outprefix} \
            --num_generations ${num_generations} \
            --migration_rate ${migration_rate} \
            --num_cells_downsample ${num_cells_downsample} \
            --max_anatomical_sites ${max_anatomical_sites} \
            --seed ${seed}

        done
    done
done
