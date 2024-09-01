#!/bin/bash

treefile=$1
#figtree_template=$2
### Below works with script called from home directory to avoid having to always input figtree template path
figtree_template="./inputs/template_figtree_block.tree"

#figtree="/home/staklins/bin/FigTree_v1.4.4/lib/figtree.jar"
figtree="$(which figtree.jar)"

# combine figtree template with tree file
figtree_treefile="${treefile/.tree/_figtree.tree}"
outfile="${treefile/.tree/_figtree.pdf}"
cat $treefile $figtree_template > $figtree_treefile

java -jar $figtree -graphic PDF $figtree_treefile $outfile
