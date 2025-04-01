#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import pickle

# meant to be used for combining multiple precision/recall pickle files into one, such as with the variable rates simulated data where data was saved per rate pair
# files = sys.argv[1]
# outfile = sys.argv[2]

# testing
files = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig7_mut005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig4_mut0005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig4_mut0025/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig5_mut0005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig4_mut001/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig6_mut001/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig4_mut01/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig6_mut01/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig7_mut001/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig5_mut001/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig5_mut01/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig5_mut005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig6_mut0025/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig7_mut0025/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig6_mut0005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig6_mut005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig4_mut005/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig5_mut0025/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig7_mut01/precision_recall_vars.pkl,/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/mig7_mut0005/precision_recall_vars.pkl"
outfile = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/variable_migration_and_mutation_rates_repeat_origin_scaling_implemented_10_25_24_data_from_8_19_24/precision_recall_variable_rates/precision_recall_vars_all.pkl"

fs = files.split(",")

# Initialize overall variables
overall_machina_precisions = []
overall_machina_recalls = []
overall_metient_precisions = []
overall_metient_recalls = []
overall_random_precisions = []
overall_random_recalls = []
overall_consensus_precisions = []
overall_consensus_recalls = []
overall_parsimony_precisions = []
overall_parsimony_recalls = []
overall_all_thresh_df = pd.DataFrame()
# overall_pathfinder_all_thresh_df = pd.DataFrame()

# Process each file and append data to overall variables
for f in fs:
    with open(f, "rb") as file:
        (
            machina_precisions,
            machina_recalls,
            metient_precisions,
            metient_recalls,
            random_precisions,
            random_recalls,
            consensus_precisions,
            consensus_recalls,
            parsimony_precisions,
            parsimony_recalls,
            all_thresh_df,
        ) = pickle.load(file)
        # machina_precisions, machina_recalls, metient_precisions, metient_recalls, random_precisions, random_recalls, consensus_precisions, consensus_recalls, parsimony_precisions, parsimony_recalls, all_thresh_df, pathfinder_all_thresh_df = pickle.load(file)

        overall_machina_precisions.extend(machina_precisions)
        overall_machina_recalls.extend(machina_recalls)
        overall_metient_precisions.extend(metient_precisions)
        overall_metient_recalls.extend(metient_recalls)
        overall_random_precisions.extend(random_precisions)
        overall_random_recalls.extend(random_recalls)
        overall_consensus_precisions.extend(consensus_precisions)
        overall_consensus_recalls.extend(consensus_recalls)
        overall_parsimony_precisions.extend(parsimony_precisions)
        overall_parsimony_recalls.extend(parsimony_recalls)
        overall_all_thresh_df = pd.concat(
            [overall_all_thresh_df, all_thresh_df], ignore_index=True
        )
        # overall_pathfinder_all_thresh_df = pd.concat([overall_pathfinder_all_thresh_df, pathfinder_all_thresh_df], ignore_index=True)

# Save the overall variables to the outfile
with open(outfile, "wb") as out_file:
    # pickle.dump((overall_machina_precisions, overall_machina_recalls, overall_metient_precisions, overall_metient_recalls, overall_random_precisions, overall_random_recalls, overall_consensus_precisions, overall_consensus_recalls, overall_parsimony_precisions, overall_parsimony_recalls, overall_all_thresh_df, overall_pathfinder_all_thresh_df), out_file)
    pickle.dump(
        (
            overall_machina_precisions,
            overall_machina_recalls,
            overall_metient_precisions,
            overall_metient_recalls,
            overall_random_precisions,
            overall_random_recalls,
            overall_consensus_precisions,
            overall_consensus_recalls,
            overall_parsimony_precisions,
            overall_parsimony_recalls,
            overall_all_thresh_df,
        ),
        out_file,
    )
