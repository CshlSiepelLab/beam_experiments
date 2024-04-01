#!/usr/bin/env python3
import os, sys
import pandas as pd
import cassiopeia as cas


def main():
    character_matrix_tsv = sys.argv[1]
    character_matrix_tsv = "/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/sim_data_barcodes_modifiedTTPmachina_3_29_24/mS/987/987_indel_character_matrix.tsv"

    # get dir path to input file for output to the same dir
    outprefix = os.path.dirname(character_matrix_tsv)

    # read in final matrix
    final_matrix = pd.read_csv(character_matrix_tsv, sep = '\t', index_col=0)

    # solve cassiopeia greedy
    reconstructed_tree = cas.data.CassiopeiaTree(character_matrix = final_matrix, missing_state_indicator = -1)
    greedy_solver = cas.solver.VanillaGreedySolver()
    greedy_solver.solve(reconstructed_tree)

    out_tree_infer = outprefix + "/cassiopeia_greedy_inferred.nwk"
    with open(out_tree_infer,'w') as it:
        it.write(reconstructed_tree.get_newick())

if __name__ == "__main__":
    main()