#!/bin/bash

indel_matrix_file=$1
tissues_tsv_file=$2
generations=$3
outdir=$4

# assumes the output dir is named with the sim number
outname=$(echo $outdir | awk -F'/' '{print $NF}')

# write tip traits to new csv
sed 's/ /,/g' $tissues_tsv_file > ${outdir}/${outname}_tip_tissues.csv

# write date trait
sed 's/ /,/g' $tissues_tsv_file | cut -d',' -f1 | paste -d',' - <(yes $generations | head -n $(wc -l < $tissues_tsv_file)) > ${outdir}/${outname}_date_traits.csv

# write fasta for tips based on input indel matrix
all_seqs=""
while IFS=$'\t' read -r -a row; do
    seq_name="${row[0]}"
    sequence="${row[@]:1}"
    sequence_csv=$(echo $sequence | sed 's/ /,/g' | sed 's/-1/0/g')
    all_seqs+=">$seq_name\n$sequence_csv\n"
done < <(tail -n +2 "$indel_matrix_file")
echo -e $all_seqs > ${outdir}/${outname}.fasta

input_names[$dataset]+="${outname},"
