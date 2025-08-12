
import sys
from ete3 import Tree

newick = sys.argv[1]
id = sys.argv[2]
outdir = sys.argv[3]

# Read the BEAST tree with tissue annotations
tree = Tree(newick, format=1)

# convert the annotation to name labels
tree_beast_annotation = tree.copy()
i = 0
for node in tree_beast_annotation.traverse():
    name = node.name
    if name.startswith("[&location="):
        name = f"node{i}{name}"
        i += 1
    name = name.replace('[&location="', "_")
    name = name.replace('"]', "")
    node.name = name

# remove annotations
tree_no_annotations = tree_beast_annotation.copy()
for node in tree_no_annotations.traverse():
    node.name = node.name.split("_")[0]

# get a tsv file of the annotations for the tips
with open(f"{outdir}/{id}_tip_tissues.tsv", "w") as f:
    for node in tree_beast_annotation.iter_leaves():
        name, tissue = node.name.split("_")
        f.write(f"{name}\t{tissue}\n")


# Write the trees to files
tree_no_annotations.write(outfile=f"{outdir}/{id}.newick", format=8)
tree_beast_annotation.write(outfile=f"{outdir}/{id}_beast.newick", format=8)
