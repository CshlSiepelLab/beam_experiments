#!/bin/bash

### This script interfaces with the barcode_simulator repo
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
# conda activate simulate

MIGRATION_MATRIX="NA"

if [[ $# -eq 0 ]] ; then
    echo "Usage: sim_wrapper.sh --design <BC10v0 or TAPE or RANDOM> --sites <num_sites> --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> [--migrationrate <float>] [--migration <filepath>] [--sites <int>]"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d|--design) DESIGN="$2"; shift ;;
        -o|--out) NAME="$2"; shift ;;
        -n|--sites) NUM_SITES="$2"; shift ;;
        -m|--mutrate) MUTRATE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;
        -mr|--migrationrate) MIGRATION_RATE="$2"; shift;;
        -mm|--migration) MIGRATION_MATRIX="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_wrapper.sh --design <BC10v0 or TAPE> --sites <num_sites> --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> [--migrationrate <float>] [--migration <filepath>] [--sites <int>]"; exit 1 ;;
    esac
    shift
done

if [ ! -f "../barcode_simulator/scripts/simulator/simulator.py" ]
then
    echo "Script ../barcode_simulator/scripts/simulator/simulator.py not found. Exiting!"
    exit
fi

outputdir="sim_results_${NAME}/"
mkdir ${outputdir}
outprefix="${outputdir}/${NAME}"

../barcode_simulator/scripts/simulator/simulator.py ${outprefix} ${DESIGN} ${NUM_SITES} ${MUTRATE} ${NUM_SAMPLES} ${MIGRATION_RATE} ${MIGRATION_MATRIX}

# remove fasta files since methods will work directly from mutation matrix, assumign this can be deduced from sequencing data accurately by other methods
rm ${outputdir}/*.fa
