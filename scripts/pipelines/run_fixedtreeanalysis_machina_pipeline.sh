#!/bin/bash

# This script will simulate true trees with groupd truth tissue location data, and then run both BEAST FixedTreeAnalysis and MACHINA to then compare the results of accuracy of internal node tissue location predictions and runtime

conda activate simulate

# Specify the number of simulations to be run for ground truth trees with migration data
num_trees=10

for ((i=1; i<=${num_trees}; i++))
do
./scripts/sim_wrapper.sh --design RANDOM --out sim${i} --sites 10 --mutrate 1.0 --samples 50 --migration inputs/test_migration_prob_matrix.csv
done
