#!/bin/bash

# inputdir=$1
inputdir="nested_sampling_uniform_vs_non_uniform_pR_rates_6_12_24/non_uniform_pR_6_6_24"

outputdir="${inputdir}/proper_joint_beast_inference"
mkdir $outputdir

files=$(find $inputdir -type f -name *_indel_character_matrix.tsv)

# total time from root to tips
num_generations=250

########################
### Format input files
########################
all_names=""
for indel_matrix_file in $files; do
    dirname=$(dirname $indel_matrix_file)
    tissues_tsv_file=${dirname}/cell_tree_*[0-9].labeling
    outname1=$(basename $dirname)
    dirname2=$(dirname $dirname)
    outname2=$(basename $dirname2)
    outname=${outname2}_${outname1}

    # write tip traits to new csv
    sed 's/ /,/g' $tissues_tsv_file > ${outputdir}/${outname}_tip_tissues.csv

    # write date trait
    sed 's/ /,/g' $tissues_tsv_file | cut -d',' -f1 | paste -d',' - <(yes $num_generations | head -n $(wc -l < $tissues_tsv_file)) > ${outputdir}/${outname}_date_traits.csv

    # write fasta for tips based on input indel matrix
    all_seqs=""
    while IFS=$'\t' read -r -a row; do
        seq_name="${row[0]}"
        sequence="${row[@]:1}"
        sequence_csv=$(echo $sequence | sed 's/ /,/g' | sed 's/-1/0/g')
        all_seqs+=">$seq_name\n$sequence_csv\n"
    done < <(tail -n +2 "$indel_matrix_file")

    echo -e $all_seqs > ${outputdir}/${outname}.fasta
    all_names+="${outname},"
done
all_names=${all_names%,}

################################################
### Run BEAST joint inference across all files
################################################
# specify tissue CTMC model
metastabayes_jar="/grid/siepel/home_norepl/staklins/metastabayes/metastabayes.jar"
xml="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/proper_joint_inference_beast.xml"
models=(oneRate threeRates sym asym)

for model in ${models[@]}; do

if [[ "$model" == "sym" ]]; then
    symmetric="true"
    spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel"
    num_rates="190"
elif [[ "$model" == "asym" ]]; then
    symmetric="false"
    spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel"
    num_rates="380"
elif [[ "$model" == "oneRate" ]]; then
    symmetric="true"
    spec="metastabayes.substitutionmodel.OneRateAllTissues"
    num_rates="1"
elif [[ "$model" == "threeRates" ]]; then
    symmetric="false"
    spec="metastabayes.substitutionmodel.ThreeRatesForSeedingRoutes"
    num_rates="3"
fi

model_dir="${outputdir}/${model}"
mkdir $model_dir
model_xml="${model_dir}/proper_joint_inference_beast_nested_sampling.xml"
cp $xml $model_xml

# replace necessary portions of the xml to run nested sampling
existing_mcmc='<run id="mcmc" spec="MCMC"'
replace_mcmc='<run id="mcmc" spec="beast.gss.NS" chainLength="1000000" particleCount="1" subChainLength="10000" epsilon="1e-10">'
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$model_xml"

existing_logger='<logger id="tracelog" spec="Logger'
replace_logger='<logger id="tracelog" spec="NSLogger'
sed -i "s|$existing_logger|$replace_logger|" "$model_xml"

beast_log="${model_dir}/joint_inference_beast_terminal_time.log"
cmd_file="${model_dir}/run_${model}.sh"

echo -e "java -Xmx5g -jar ${metastabayes_jar} -overwrite -working \
-D "inputNames=${all_names}" \
-D "generations=${num_generations}" \
-D "fileDir=${outputdir}" \
-D "traitModelSpec=${spec}" \
-D "symmetric=${symmetric}" \
-D "numRates=${num_rates}" \
$model_xml > $beast_log" > $cmd_file


qsub -cwd -l m_mem_free=5G -o ${cmd_file}.out -e ${cmd_file}.err $cmd_file

done




