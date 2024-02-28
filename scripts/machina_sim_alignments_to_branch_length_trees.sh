#!/bin/bash

### This is a basic bash script to convert MACHINA sim data alignments into branch length resolved trees using mega-cc similar to PathFinder's approach

files=$(find machina_data/sims/ -type f -name reads_*.fas)

for file in $files;
do

outputdir=$(dirname "$file")

python scripts/pathfinder/pathfinder.py $file -o $outputdir

#rm -r ${outputdir}/M11CC_Out
mv ${outputdir}/scratch* ${outputdir}/megacc

done
