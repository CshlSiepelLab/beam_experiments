
import sys
from beam_sup.tree_inference import infer_parsimony_tree_cassiopeia_greedy

character_matrix_tsv = sys.argv[1]
outprefix = sys.argv[2]

infer_parsimony_tree_cassiopeia_greedy(character_matrix_tsv, outprefix)

