#!/usr/bin/env python

import sys
import ast
import os
import ete3
from Bio import Phylo
from io import StringIO


machina = sys.argv[1]
metient = sys.argv[2]
beam = sys.argv[3]
parsimony = sys.argv[4]
primaryTissue = sys.argv[5]
outdir = sys.argv[6]
name = sys.argv[7] 

# # testing
# machina = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/machina/MMUS1469/CP10/PRL-G-PRL-R.tree"
# metient = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/metient/MMUS1469/CP10/MMUS1469_CP10_PRL_migration_graphs.txt"
# beam = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/beam/MMUS1469/CP10/combined.trees"
# parsimony = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/serio_prostate_cancer_data_11_18_24_asv_cutoff_50/parsimony_tissue_inference/MMUS1469/CP10/all_parsimony_solutions.txt"
# primaryTissue = "PRL"
# outdir = "./"
# name = "test"


def getMigrationComigrationCountsFromNewick(newick):
    tree = ete3.Tree(newick, format=8)
    migration_count = 0
    comigration_count = 0
    met_to_met = False
    reseeding = False
    edges = []
    comigration_edges_already_checked = []

    for node in tree.traverse():
        if node.is_root():
            continue
        if node.up.is_root():
            parent_tissue = primaryTissue
        else:
            parent_tissue = node.up.name.split("_")[-1]
        child_tissue = node.name.split("_")[-1]
        if parent_tissue != child_tissue:
            migration_count += 1
            edge = f"{parent_tissue}_{child_tissue}"
            if edge in edges and edge not in comigration_edges_already_checked:
                comigration_count += 1
                comigration_edges_already_checked.append(edge)
            if not met_to_met and parent_tissue != primaryTissue and child_tissue != primaryTissue:
                met_to_met = True
            if not reseeding and child_tissue == primaryTissue:
                reseeding = True
            edges.append(edge)
    return migration_count, comigration_count, met_to_met, reseeding

# get machina counts and boolean topology classifications
if os.path.isfile(machina):
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
            if not machina_met_to_met and source != primaryTissue and recipient != primaryTissue:
                machina_met_to_met = True
            if not machina_reseeding and recipient == primaryTissue:
                machina_reseeding = True

    machina_migration_count = sum([value for key, value in machina_migrations.items()])
    machina_comigration_count = len(machina_migrations.keys())
else:
    machina_migration_count = float('nan')
    machina_comigration_count = float('nan')
    machina_met_to_met = float('nan')
    machina_reseeding = float('nan')

# get metient counts and boolean topology classifications
if os.path.isfile(metient):
    metient_loss_graphs = {}
    with open(metient, "r") as f:
        next(f)  # Skip the first line
        i=1
        for line in f.readlines():
            loss, graph = line.strip().split("\t")
            metient_loss_graphs[i] = graph
            i+=1

    all_migration_counts = []
    all_comigration_counts = []
    met_to_met_counts = []
    reseeding_counts = []

    for loss, graph in metient_loss_graphs.items():
        graph_dict = ast.literal_eval(graph)

        migration_count = 0
        comigration_count = 0
        met_to_met = False
        reseeding = False

        for key, value in graph_dict.items():
            for k, v in value.items():
                if v == 0:
                    continue
                else:
                    migration_count += int(v)
                    comigration_count += 1
                    if not met_to_met and key != primaryTissue and k != primaryTissue:
                        met_to_met = True
                    if not reseeding and k == primaryTissue:
                        reseeding = True

        all_migration_counts.append(migration_count)
        all_comigration_counts.append(comigration_count)
        met_to_met_counts.append(met_to_met)
        reseeding_counts.append(reseeding)

    metient_migration_count = sum(all_migration_counts) / len(all_migration_counts)
    metient_comigration_count = sum(all_comigration_counts) / len(all_comigration_counts)
    metient_met_to_met = max(set(met_to_met_counts), key=met_to_met_counts.count)
    metient_reseeding = max(set(reseeding_counts), key=reseeding_counts.count)
else:
    metient_migration_count = float('nan')
    metient_comigration_count = float('nan')
    metient_met_to_met = float('nan')
    metient_reseeding = float('nan')

# get parsimony counts and boolean topology classifications
if os.path.isfile(parsimony):
    parsimony_solutions = []
    with open(parsimony, "r") as f:
        for line in f.readlines():
            parsimony_solutions.append(line.strip())
    
    migration_counts = []
    comigration_counts = []
    met_to_mets = []
    reseedings = []
    for newick in parsimony_solutions:
        mig_count, comig_count, met_met, reseeding = getMigrationComigrationCountsFromNewick(newick)
        migration_counts.append(mig_count)
        comigration_counts.append(comig_count)
        met_to_mets.append(met_met)
        reseedings.append(reseeding)
    
    parsimony_migration_count = sum(migration_counts) / len(migration_counts)
    parsimony_comigration_count = sum(comigration_counts) / len(comigration_counts)
    parsimony_met_to_met = max(set(met_to_mets), key=met_to_mets.count)
    parsimony_reseeding = max(set(reseedings), key=reseedings.count)

# get beam counts and boolean topology classifications
if os.path.isfile(beam):
    beam_migration_counts = []
    beam_comigration_counts = []
    beam_met_to_met_counts = []
    beam_reseeding_counts = []

    trees = Phylo.parse(beam, 'nexus')

    num_trees = len([str(tree) for tree in trees])
    burnin = int(num_trees * 0.1)

    i = 0
    j = 0
    # Need to remake the generator object since it disappears once it is called above to calculate the burnin amount
    trees = Phylo.parse(beam, 'nexus')
    for tree in trees:
        i += 1

        # discard burnin
        if i < burnin:
            continue

        tree.rooted = True  # ensure the tree is rooted
        tree.format = 'newick'  # set the format to newick
        for clade in tree.find_clades():
            if hasattr(clade, 'comment'):
                location = clade.comment.split('location="')[1].split('"')[0]
            if clade.name is None:
                clade.name = f"node{j}_{location}"
                j += 1
            else:
                clade.name = f"{clade.name}_{location}"
            if hasattr(clade, 'comment'):
                del clade.comment
            if hasattr(clade, 'branch_length') and clade.branch_length is not None:
                del clade.branch_length
        
        newick_io = StringIO()
        Phylo.write(tree, newick_io, 'newick', plain=True)
        newick_io.seek(0)
        bnewick = newick_io.read()

        mig_count, comig_count, met_met, reseeding = getMigrationComigrationCountsFromNewick(bnewick)
        beam_migration_counts.append(mig_count)
        beam_comigration_counts.append(comig_count)
        beam_met_to_met_counts.append(met_met)
        beam_reseeding_counts.append(reseeding)


    beam_migration_count = sum(beam_migration_counts) / len(beam_migration_counts)
    beam_comigration_count = sum(beam_comigration_counts) / len(beam_comigration_counts)
    beam_met_to_met = max(set(beam_met_to_met_counts), key=beam_met_to_met_counts.count)
    beam_reseeding = max(set(beam_reseeding_counts), key=beam_reseeding_counts.count)

# write migration and co-migration counts to file (Note: not temporally consistent co-migrations since the counts are taken from the graphs alone)
with open(f"{outdir}/{name}_migration_counts.csv", "w") as f:
    f.write("name,"
            "machina_migrations,machina_comigrations,"
            "metient_migrations,metient_comigrations,"
            "parsimony_migrations,parsimony_comigrations,"
            "beam_migrations,beam_comigrations\n"
            f"{name},"
            f"{machina_migration_count},{machina_comigration_count},"
            f"{metient_migration_count},{metient_comigration_count},"
            f"{parsimony_migration_count},{parsimony_comigration_count},"
            f"{beam_migration_count},{beam_comigration_count}")

with open(f"{outdir}/{name}_topology_classifications.csv", "w") as f:
    f.write("name,"
            "machina_met_to_met,machina_reseeding,"
            "metient_met_to_met,metient_reseeding,"
            "parsimony_met_to_met,parsimony_reseeding,"
            "beam_met_to_met,beam_reseeding\n"
            f"{name},"
            f"{machina_met_to_met},{machina_reseeding},"
            f"{metient_met_to_met},{metient_reseeding},"
            f"{parsimony_met_to_met},{parsimony_reseeding},"
            f"{beam_met_to_met},{beam_reseeding}")

