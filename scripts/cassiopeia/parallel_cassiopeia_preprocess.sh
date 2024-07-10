#!/bin/bash

outdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/preprocess_cassiopeia_7_9_24"
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
    qsub -cwd -l m_mem_free=2G -pe threads 10 -o $outdir/logs/$sra_file_id.log -e $outdir/logs/$sra_file_id.err $outdir/cmds/$sra_file_id.sh
done



############################################################################################################
# used below to resubmit failed jobs with more memory
############################################################################################################

# to get IDs without the final step run due to seg fault
# files=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/preprocess_cassiopeia_6_26_24/logs/ -type f -name *.err)
# errors=()
# for file in $files; do
#     echo $file
#     id=$(basename $file | cut -d "." -f 1)
#     if grep -q "Plotting filtered lineage group pivot table heatmap" $file; then
#         echo $id "completed"
#     else
#         errors+=($id)
#     fi
# done

# # or to use all IDs
# files=$(find /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/yang_2022_real_data/preprocess_cassiopeia_6_26_24/cmds/ -type f -name *.sh)
# errors=()
# for file in $files; do
#     echo $file
#     id=$(basename $file | cut -d "." -f 1)
#     errors+=($id)
# done

# # or to use only the filtered mice ids
# errors=(SRR17885790 SRR17885791 SRR17885792 SRR17885793 SRR17885797 SRR17885799 SRR17885819 SRR17885820 SRR17885822 SRR17885823 SRR17885824 SRR17885825 SRR17885828 SRR17885829 SRR17885834 SRR17885835 SRR17885839)

# # to rerun error jobs
# for error in ${errors[@]}; do
#         rm -r $outdir/logs/$error*
#         rm -r $outdir/$error
#         qsub -cwd -l m_mem_free=20G -pe threads 25 -o $outdir/logs/$error.log -e $outdir/logs/$error.err $outdir/cmds/$error.sh
# done

