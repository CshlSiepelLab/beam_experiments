#!/bin/bash

# This script sets up the input/output files for the nested sampling run through feast io for beast xml formatting across datasets

newick_file=$1
tissues_tsv_file=$2

# # # testing
# newick_file="machina_sims_proper_nested_sampling_5_29_24/m8/mS/seed5.nwk"
# tissues_tsv_file="machina_sims_proper_nested_sampling_5_29_24/m8/mS/seed5.labeling"

# get all tip names from the newick file by pattern matching for newick with only tip names and all branch lengths
tip_names=$(grep -o '[,(][^,:]*:' $newick_file | tr -d ',:(')

# use those names
tip_traits=""
for tip in $tip_names; do
tip_traits+="$(grep $tip $tissues_tsv_file | tr -d ' ' | sed 's/\t/,/g')\n"
done

# write tip traits to new csv
outfile_traits=${tissues_tsv_file//.labeling/_tips_only.csv}
echo -e $tip_traits > $outfile_traits

# write fake fasta for tip names
outfile_fasta=${tissues_tsv_file//.labeling/.fasta}
for name in $tip_names; do
echo -e ">$name\n?" >> $outfile_fasta
done

# reformat newick file to remove extra bracket and false root branch length *****BASED ON MACHINA SIM DATA ONLY*****
outfile_newick=${newick_file//.nwk/_reformatted.nwk}
cat $newick_file | sed '1s/^(//' $newick_file | sed 's/):[0-9]*;$/;/' > $outfile_newick
