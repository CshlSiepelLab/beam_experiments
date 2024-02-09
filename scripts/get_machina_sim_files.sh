#!/bin/bash

machina_sims="/local/storage/no-backup/staklins-scratch/machina/data/sims/m8"
destination="/local/storage/no-backup/staklins-scratch/bayesian_phylogenetic_metastasis/machina_m8_sim_data"

find $machina_sims -type f -name 'T_*.tree' -exec cp {} $destination \;
find $machina_sims -type f -name 'T_*.vertex.labeling' -exec cp {} $destination \;
