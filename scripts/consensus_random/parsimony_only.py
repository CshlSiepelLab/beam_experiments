
import sys
import pandas as pd
from ete3 import Tree

from beam_sup.migration_inference import label_tissues_fitch_parsimony


# User inputs
tree_file = sys.argv[1]  # newick file
leaf_tissues_tsv = sys.argv[
    2
]  # tsv file with tip cell names and tissue labels as columns, and no header
outdir = sys.argv[3]  # where to write the output
primary_tissue = sys.argv[4]  # the known tissue label of the root node
threshold_num_solutions = int(
    sys.argv[5]
)  # the maximum number of possible solutions to enumerate them all in the output


tree = Tree(tree_file, format=8)
tissues_df = pd.read_csv(
    leaf_tissues_tsv,
    sep=r"\s+",
    header=None,
    names=["cell", "tissue"],
    dtype={"cell": str, "tissue": str},
)

# Force the tree to be binary (ie. resolve polytomies)
tree.resolve_polytomy()

# get results
random_parsimony_tree, all_solutions = label_tissues_fitch_parsimony(
    tree, 
    tissues_df, 
    primary_tissue,
    threshold_num_solutions
)

parsimony_output = outdir + "/parsimony_tissues_random.nwk"
random_parsimony_tree.write(outfile=parsimony_output, format=8, format_root_node=True)

if len(all_solutions) > 0:
    # Write all solutions to a file
    for i, solution in enumerate(all_solutions):
        solution_output = f"{outdir}/parsimony_tissues_all_solutions_{i+1}.nwk"
        solution.write(outfile=solution_output, format=8, format_root_node=True)
