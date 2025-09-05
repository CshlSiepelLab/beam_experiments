
import sys

from beam_sup.tree_utils import get_mig_comig_counts_and_topologies_from_nwk

newick_file = sys.argv[1]
primary_tissue = sys.argv[2]
output_file = sys.argv[3]

migration_count, comigration_count, num_multiedges, met_to_met, primary_reseeding, clonality = get_mig_comig_counts_and_topologies_from_nwk(
    newick_file,
    primary_tissue
)

with open(output_file, "w") as f:
    f.write(f"met_to_met,primary_reseeding\n")
    f.write(f"{met_to_met},{primary_reseeding}\n")
