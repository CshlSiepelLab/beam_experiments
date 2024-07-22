#!/bin/bash

# specify the path to the directory where the results will be stored
outdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/results/model_selection_results_7_11_24"

# specify the path to the directory where the data is
dataset_dirs="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/model_selection_data_7_8_24"
dataset_names=$(find $dataset_dirs -mindepth 1 -maxdepth 1 -type d | awk -F'/' '{print $NF}')


# make map  to track sim names for each dataset
unset input_names
declare -A input_names


generations=250

for dataset in $dataset_names; do
    echo $dataset

    # set empty input names
    input_names[$dataset]=""

    dataset_outdir=$outdir/$dataset
    mkdir -p $dataset_outdir

    # Format input files for feast read in to beast proper joint inference
    indel_matrix_files=$(find $dataset_dirs/$dataset -type f -name *_indel_character_matrix.tsv)

    for indel_matrix_file in $indel_matrix_files; do
        dirname=$(dirname $indel_matrix_file)
        outname=$(echo $dirname | awk -F'/' '{print $NF}')
        tissues_tsv_file=${dirname}/cell_tree_*[0-9].labeling

        # write tip traits to new csv
        sed 's/ /,/g' $tissues_tsv_file > ${dataset_outdir}/${outname}_tip_tissues.csv

        # write date trait
        sed 's/ /,/g' $tissues_tsv_file | cut -d',' -f1 | paste -d',' - <(yes $generations | head -n $(wc -l < $tissues_tsv_file)) > ${dataset_outdir}/${outname}_date_traits.csv

        # write fasta for tips based on input indel matrix
        all_seqs=""
        while IFS=$'\t' read -r -a row; do
            seq_name="${row[0]}"
            sequence="${row[@]:1}"
            sequence_csv=$(echo $sequence | sed 's/ /,/g' | sed 's/-1/0/g')
            all_seqs+=">$seq_name\n$sequence_csv\n"
        done < <(tail -n +2 "$indel_matrix_file")
        echo -e $all_seqs > ${dataset_outdir}/${outname}.fasta

        input_names[$dataset]+="${outname},"
    done
    # remove trailing comma from input_names
    input_names[$dataset]=${input_names[$dataset]%,}
done



# setup directories for xml files and run BEAST2 nested sampling proper joint inference
models=(model1 model2 model3 model4 model5)
template_xml="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/proper_nested_sampling_template.xml"
metastabayes_jar="/grid/siepel/home_norepl/staklins/metastabayes/metastabayes.jar"

for model1 in ${models[@]}; do
    sim_names=${input_names[$model1]}
    file_dir=$outdir/$model1
    for model2 in ${models[@]}; do
        dir="$outdir/${model1}/${model2}_ns_results"
        mkdir -p $dir
        cp $template_xml $dir
        matrix_structure=$(cat /grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/model_selection_transition_matrices/${model2}_inference_matrixStructure.csv | tr '\n' ' ')

        # assumes that the unique rates are numbered starting from 0 (not inclusing the use of 0 for the irrelevant diagonal rates which are to be replaced automatically with metastabayes java code)
        num_rates=$(echo "$matrix_structure" | tr ' ' '\n' | sort | uniq | wc -l)

        run_file=$dir/run_${model1}_data_${model2}_inference.sh
        echo -e "java -jar $metastabayes_jar -threads 17 -overwrite -working -D inputNames=$sim_names -D fileDir=$file_dir -D matrixStructure='$matrix_structure' -D numRates=$num_rates $dir/proper_nested_sampling_template.xml > $dir/run_${model1}_data_${model2}_inference_terminal.log" > $dir/run_${model1}_data_${model2}_inference.sh
        chmod u+x $run_file
        qsub -cwd -l m_mem_free=5G -pe threads 17 -e $dir/cluster.log -o $dir/cluster.log $run_file
    done
done
