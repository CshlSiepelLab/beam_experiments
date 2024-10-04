#!/usr/bin/env python

import sys
import ast

machina = sys.argv[1]
metient = sys.argv[2]
metastabayes = sys.argv[3]
primaryTissue = sys.argv[4]
outdir = sys.argv[5]
name = sys.argv[6] 

# # testing
# machina = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/machina/MMUS1457/CP02/PRL-G-PRL-R.tree"
# metient = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/metient/MMUS1457/CP02/MMUS1457_CP02_PRL_migration_graphs.txt"
# metastabayes = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/metastabayes/MMUS1457/CP02/posterior_prob_graph.csv"
# primaryTissue = "PRL"
# outdir = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/asv50_ryan_prostate_cancer_data_9_5_24/compare_migration_counts"
# name = "test"

# set threshold for metastabayes
prob_consensus_threshold = 0.7

# get machina counts and boolean topology classifications
machina_met_to_met = False
machina_reseeding = False
machina_migrations = {}
with open(machina, "r") as f:
    for line in f.readlines():
        migration = "_".join(line.strip().split(" "))
        if migration in machina_migrations:
            machina_migrations[migration] += 1
        else:
            machina_migrations[migration] = 1
        source, recipient = migration.split("_")
        if not machina_met_to_met and source != primaryTissue:
            machina_met_to_met = True
        if not machina_reseeding and recipient == primaryTissue:
            machina_reseeding = True

machina_migration_count = sum([value for key, value in machina_migrations.items()])
machina_comigration_count = len(machina_migrations.keys())

# get metient counts and boolean topology classifications
metient_loss_graphs = {}
with open(metient, "r") as f:
    next(f)  # Skip the first line
    for line in f.readlines():
        loss, graph = line.strip().split("\t")
        metient_loss_graphs[loss] = graph

lowest_loss_graph_dict = ast.literal_eval(metient_loss_graphs[min(metient_loss_graphs.keys(), key=float)])

metient_migration_count = 0
metient_comigration_count = 0
metient_met_to_met = False
metient_reseeding = False
for key, value in lowest_loss_graph_dict.items():
    for k, v in value.items():
        if v == 0:
            continue
        else:
            metient_migration_count += int(v)
            metient_comigration_count += 1
            if not metient_met_to_met and key != primaryTissue:
                metient_met_to_met = True
            if not metient_reseeding and k == primaryTissue:
                metient_reseeding = True


# get metastabayes counts and boolean topology classifications
metastabayes_migrations = []
with open(metastabayes, "r") as f:
    for line in f.readlines():
        migration, probability = line.strip().split(",")
        if float(probability) >= prob_consensus_threshold:
            metastabayes_migrations.append("_".join(migration.split("_")[0:2]))

metastabayes_met_to_met = False
metastabayes_reseeding = False
done = []
metastabayes_migration_count = 0
metastabayes_comigration_count = 0
for migration in metastabayes_migrations:
    if migration not in done:
        metastabayes_comigration_count += 1
        done.append(migration)

        s, r = migration.split("_")
        if not metastabayes_met_to_met and s != primaryTissue:
                metastabayes_met_to_met = True
        if not metastabayes_reseeding and r == primaryTissue:
                metastabayes_reseeding = True
    metastabayes_migration_count += 1

# write migration and co-migration counts to file (Note: not temporally consistent co-migrations since the counts are taken from the graphs alone)
with open(f"{outdir}/{name}_migration_counts.csv", "w") as f:
    f.write("name,"
            "machina_migrations,machina_comigrations,"
            "metient_migrations,metient_comigrations,"
            "metastabayes_migrations,metastabayes_comigrations\n"
            f"{name},"
            f"{machina_migration_count},{machina_comigration_count},"
            f"{metient_migration_count},{metient_comigration_count},"
            f"{metastabayes_migration_count},{metastabayes_comigration_count}")

with open(f"{outdir}/{name}_topology_classifications.csv", "w") as f:
    f.write("name,"
            "machina_met_to_met,machina_reseeding,"
            "metient_met_to_met,metient_reseeding,"
            "metastabayes_met_to_met,metastabayes_reseeding\n"
            f"{name},"
            f"{machina_met_to_met},{machina_reseeding},"
            f"{metient_met_to_met},{metient_reseeding},"
            f"{metastabayes_met_to_met},{metastabayes_reseeding}")

