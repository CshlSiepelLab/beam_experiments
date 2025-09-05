
import sys
import ast
import os
import ete3
from Bio import Phylo
from io import StringIO

from beam_sup.tree_utils import get_mig_comig_counts_and_topologies_from_nwk


sim_name = sys.argv[1]
newick = sys.argv[2]
primary_tissue = sys.argv[3]
outfile = sys.argv[4]

migration_count, comigration_count, num_multiedges, met_to_met, reseeding, clonality = (
    get_mig_comig_counts_and_topologies_from_nwk(newick, primary_tissue)
)

with open(outfile, "a") as f:
    f.write(
        f"{sim_name},{migration_count},{comigration_count},{num_multiedges},{met_to_met},{reseeding},{clonality}\n"
    )
