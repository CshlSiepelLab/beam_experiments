#!/bin/bash

inputdir=$1
# inputdir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/data/new_simulator_unifromTransitionProbs_6_6_24"

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
model="oneRate"

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

metastabayes_jar="/grid/siepel/home_norepl/staklins/metastabayes/metastabayes.jar"
xml="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/proper_joint_inference_beast.xml"
beast_dir="${outputdir}/beast"
mkdir $beast_dir

beast_logs=()
num_chains=5
for ((i=1; i<=$num_chains; i++))
do
  beast_log="${beast_dir}/joint_inference_beast_terminal_time_${i}.log"
  beast_logs+=("$beast_log")
  iter_xml="${beast_dir}/joint_inference_beast_${i}.xml"
  cp $xml $iter_xml
  time java -Xmx5g -jar ${metastabayes_jar} -overwrite -working \
    -D "inputNames=${all_names}" \
    -D "generations=${num_generations}" \
    -D "fileDir=${outputdir}" \
    -D "traitModelSpec=${spec}" \
    -D "symmetric=${symmetric}" \
    -D "numRates=${num_rates}" \
    $iter_xml > $beast_log &
done

wait

# combine all runs for each sim trees and the overall log files
treeannotator_path=$(which treeannotator)
logcombiner_path=$(which logcombiner)

# setup stats files
stats_file="${outputdir}/proper_joint_inference_accuracy.csv"
touch $stats_file
echo "migration_topology,seed,beast_mcc_f1,beast_posterior_f1,beast_posterior_95ci_binary,true_beast_mcc_f1,true_beast_posterior_f1,true_beast_posterior_95ci_binary" > $stats_file

# combine tree files
trees_files=$(find $beast_dir -type f -name joint_inference_beast_1_tissues_*.trees)
for file in $trees_files; do
    name=$(basename $file | cut -d'_' -f6-7 | sed 's/.trees//g')
    name_trees_files=$(find $beast_dir -type f -name joint_inference_beast_*_tissues_${name}.trees)
    combined_trees="${beast_dir}/${name}_tissues_combined.trees"
    $logcombiner_path $name_trees_files -o $combined_trees
    mcc_tree=$(echo "$combined_trees" | sed 's/.trees/.tree/')
    ${treeannotator_path} -burnin 10 -topology MCC -height mean -file ${combined_trees} ${mcc_tree}

    # get migration graph performance statistics
    migration_topology=$(echo $name | cut -d'_' -f1)
    seed=$(echo $name | cut -d'_' -f2)
    true_tissue_tree="${dirname}/${migration_topology}/${seed}/*_tissue_labeled_tree.nwk"
    true_migration_graph="${dirname}/${migration_topology}/${seed}/migration_graph_seed*.csv"
    # stats for true tree
    python scripts/format_treeannotator_nexus_to_newick.py ${mcc_tree}
    beast_mcc_f1=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_tissue_tree} ${mcc_tree}.nwk | awk -F' ' '{print $3}')
    # calculate the same F1 score but for sampling all trees from the beast posterior with F1 score weighted by posterior probability
    beast_posterior=$(python scripts/migration_graph_f1_true_beast_posterior_trees.py ${true_tissue_tree} ${combined_trees})
    beast_posterior_f1=$(echo $beast_posterior | awk -F' ' '{print $3}')
    beast_posterior_95ci_binary=$(echo $beast_posterior | awk -F' ' '{print $NF}')
    # stats for true graph
    true_beast_mcc_f1=$(python scripts/migration_graph_f1_true_inferred_trees.py ${true_migration_graph} ${mcc_tree}.nwk | awk -F' ' '{print $3}')
    true_beast_posterior=$(python scripts/migration_graph_f1_true_beast_posterior_trees.py ${true_migration_graph} ${combined_trees})
    true_beast_posterior_f1=$(echo $true_beast_posterior | awk -F' ' '{print $3}')
    true_beast_posterior_95ci_binary=$(echo $true_beast_posterior | awk -F' ' '{print $NF}')
    echo -e "${migration_topology},${seed},${beast_mcc_f1},${beast_posterior_f1},${beast_posterior_95ci_binary},${true_beast_mcc_f1},${true_beast_posterior_f1},${true_beast_posterior_95ci_binary}" >> $stats_file
done

# combine log files
log_files=$(find $beast_dir -type f -name *.log)
combined_log="${beast_dir}/joint_inference_beast_combined.log"
$logcombiner_path $log_files -o $combined_log



