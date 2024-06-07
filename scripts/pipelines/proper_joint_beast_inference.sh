#!/bin/bash

inputdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/new_simulator_unifromTransitionProbs_6_6_24"

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
    # get output seed name by the two parent directories of the input file to be unique to each seed across different migration topology simulations
    dirname=$(dirname $indel_matrix_file)
    tissues_tsv_file=${dirname}/cell_tree_*[0-9].labeling
    outname1=$(basename $dirname)
    dirname2=$(dirname $dirname)
    outname2=$(basename $dirname2)
    outname=${outname2}_${outname1}

    echo $indel_matrix_file
    echo $tissues_tsv_file
    echo $outname

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

################################################
### Run BEAST joint inference across all files
################################################
# specify tissue CTMC model
model="asym"

if [[ "$model" == "sym" ]]; then
    symmetric="true"
    spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel"
    num_rates="55"
elif [[ "$model" == "asym" ]]; then
    symmetric="false"
    spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel"
    num_rates="110"
elif [[ "$model" == "oneRate" ]]; then
    symmetric="true"
    spec="metastabayes.substitutionmodel.OneRateAllTissues"
    num_rates="1"
elif [[ "$model" == "threeRates" ]]; then
    symmetric="false"
    spec="metastabayes.substitutionmodel.ThreeRatesForSeedingRoutes"
    num_rates="3"
fi

xml="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/proper_joint_inference_beast.xml"
beast_dir="${outputdir}/beast"
mkdir $beast_dir

beast_logs=()
num_chains=1
for ((i=1; i<=$num_chains; i++))
do
  beast_log="${beast_dir}/joint_inference_beast_terminal_time_${i}.log"
  beast_logs+=("$beast_log")
  iter_xml="${beast_dir}/joint_inference_beast_${i}.xml"
  cp $xml $iter_xml
  time java -Xmx5g -jar ${metastabayes_jar} -overwrite -working \
    -D "inputNames=${all_names}" \
    -D "generations=${num_generations} \
    -D "fileDir=${file_dir}" \
    -D "traitModelSpec=${spec}" \
    -D "symmetric=${symmetric}" \
    -D "numRates=${num_rates}" \
    $iter_xml > $beast_log &
done


