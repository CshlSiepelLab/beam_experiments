
import sys
import os
from ete3 import Tree
import cassiopeia as cas
import pandas as pd
import networkx as nx


def to_newick(tree: nx.DiGraph, record_branch_lengths: bool = False) -> str:
    """Converts a networkx graph to a newick string.

    Args:
        tree: A networkx tree
        record_branch_lengths: Whether to record branch lengths on the tree in
            the newick string

    Returns:
        A newick string representing the topology of the tree
    """

    def _to_newick_str(g, node):
        is_leaf = g.out_degree(node) == 0
        weight_string = ""

        if record_branch_lengths and g.in_degree(node) > 0:
            parent = list(g.predecessors(node))[0]
            weight_string = ":" + str(g[parent][node]["length"])

        _name = str(node)
        return (
            "%s" % (_name,) + weight_string
            if is_leaf
            else (
                "("
                + ",".join(_to_newick_str(g, child) for child in g.successors(node))
                + ")"
                + _name
                + weight_string
            )
        )

    root = [node for node in tree if tree.in_degree(node) == 0][0]
    return _to_newick_str(tree, root) + ";"


newick_file = sys.argv[1]
tissues_file = sys.argv[2]
outdir = sys.argv[3]

# read in newick to ete3 tree
ete_tree = Tree(newick_file, format=3)
ete_tree.name = "root"

# load the tissues to dictionary
tissues_df = pd.read_csv(
    tissues_file, header=None, index_col=0, names=["cell", "tissue"], dtype=str
)
tissues_df.index = tissues_df.index.astype(str)

# load the tree to cassiopeia object
tree = cas.data.CassiopeiaTree(tree=ete_tree, cell_meta=tissues_df)

# run fitch-hartigan to get a randomly selected parsimonious tissue labeling on the tree
fh_tree = cas.tl.fitch_hartigan(cassiopeia_tree=tree, meta_item="tissue", copy=True)

# show and save results (keep in mind that only the origin is known above the root, which the parsimony here does not consider but it should not matter)
name_map = {}
for node in fh_tree.depth_first_traverse_nodes(postorder=False):
    # attribute_names = fh_tree._CassiopeiaTree__network.nodes[node].keys()
    # print(node, list(attribute_names))
    label = fh_tree.get_attribute(node, "label")
    new_name = f"{node}_{label}"
    name_map[node] = new_name
fh_tree.relabel_nodes(name_map)

with open(f"{outdir}/cassiopeia_fitch_hartigan_result.nwk", "w") as f:
    f.write(to_newick(fh_tree._CassiopeiaTree__network))


# Optional: run fitch-count to get a transition matrix of tissue change frequencies across all parsimonious labelings
fc_matrix = cas.tl.fitch_count(cassiopeia_tree=tree, meta_item="tissue")
fc_matrix.to_csv(f"{outdir}/cassiopeia_fitch_count_result.csv")
