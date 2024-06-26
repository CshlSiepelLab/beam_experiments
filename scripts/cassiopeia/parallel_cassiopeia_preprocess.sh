#!/bin/bash

outdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/preprocess_cassiopeia_6_26_24"
mkdir -p $outdir
mkdir -p $outdir/cmds
mkdir -p $outdir/logs

unique_sra_file_ids=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/raw_fastqs -type f -name *_1.fastq.gz)

for file in $unique_sra_file_ids; do
    sra_file_id=$(basename $file | cut -d "_" -f 1)
    r1=$file
    r2=$(echo $r1 | sed 's/_1.fastq.gz/_2.fastq.gz/')
    output_dir=$outdir/$sra_file_id
    echo "python /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/scripts/cassiopeia/cassiopeia_preprocess.py $r1,$r2 $output_dir" > $outdir/cmds/$sra_file_id.sh
    qsub -cwd -l m_mem_free=2G -pe threads 5 -o $outdir/logs/$sra_file_id.log -e $outdir/logs/$sra_file_id.err $outdir/cmds/$sra_file_id.sh
done