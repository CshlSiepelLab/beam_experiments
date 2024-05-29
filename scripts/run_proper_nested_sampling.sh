#!/bin/bash

metastabayes_jar="/grid/siepel/home_norepl/staklins/metastabayes/metastabayes.jar"
# omit m8_seed394 for bug, return to fix later
# input_names="inputNames=m5_seed518,m5_seed25,m5_seed955,m5_seed0,m5_seed35,m5_seed907,m5_seed4,m5_seed865,m5_seed40,m5_seed5,m5_seed247,m5_seed1140,m5_seed31,m5_seed9,m5_seed950,m5_seed12,m5_seed76,m5_seed571,m5_seed49,m5_seed81,m5_seed62,m5_seed981,m5_seed565,m5_seed10,m5_seed2,m5_seed473,m5_seed17,m5_seed3,m5_seed538,m5_seed2155,m5_seed7,m5_seed512,m5_seed8,m5_seed32,m5_seed209,m5_seed545,m5_seed694,m5_seed534,m5_seed23,m8_seed216,m8_seed54,m8_seed0,m8_seed35,m8_seed1070,m8_seed4,m8_seed37,m8_seed5,m8_seed31,m8_seed9,m8_seed12,m8_seed76,m8_seed19,m8_seed10046,m8_seed241,m8_seed239,m8_seed69,m8_seed45,m8_seed981,m8_seed30342,m8_seed10,m8_seed2,m8_seed172,m8_seed3,m8_seed157,m8_seed905,m8_seed383,m8_seed7,m8_seed8,m8_seed243,m8_seed10157,m8_seed23"
file_dir="fileDir=/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/machina_sims_proper_nested_sampling_5_29_24"
output_dir="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/machina_sims_proper_nested_sampling_5_29_24"
input_xml="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/inputs/proper_nested_sampling_machina_sims.xml"

unset input_names
declare -A input_names
input_names["mS"]="inputNames=m5/mS/seed2,m5/mS/seed8,m5/mS/seed9,m5/mS/seed10,m5/mS/seed5,m5/mS/seed0,m5/mS/seed12,m5/mS/seed4,m5/mS/seed3,m5/mS/seed7,m8/mS/seed2,m8/mS/seed8,m8/mS/seed9,m8/mS/seed10,m8/mS/seed5,m8/mS/seed0,m8/mS/seed12,m8/mS/seed4,m8/mS/seed3,m8/mS/seed7"
input_names["S"]="inputNames=m5/S/seed62,m5/S/seed17,m5/S/seed81,m5/S/seed23,m5/S/seed49,m5/S/seed25,m5/S/seed40,m5/S/seed35,m5/S/seed31,m5/S/seed32,m8/S/seed54,m8/S/seed37,m8/S/seed2,m8/S/seed23,m8/S/seed69,m8/S/seed5,m8/S/seed0,m8/S/seed12,m8/S/seed35,m8/S/seed31"
input_names["M"]="inputNames=m5/M/seed545,m5/M/seed565,m5/M/seed512,m5/M/seed473,m5/M/seed76,m5/M/seed694,m5/M/seed209,m5/M/seed534,m5/M/seed865,m5/M/seed907,m8/M/seed239,m8/M/seed45,m8/M/seed243,m8/M/seed76,m8/M/seed172,m8/M/seed35,m8/M/seed19,m8/M/seed216,m8/M/seed7"
# omit R m8_seed394 for bug, return to fix later
input_names["R"]="inputNames=m5/R/seed538,m5/R/seed955,m5/R/seed17,m5/R/seed518,m5/R/seed571,m5/R/seed981,m5/R/seed2155,m5/R/seed247,m5/R/seed1140,m5/R/seed950,m8/R/seed30342,m8/R/seed9,m8/R/seed1070,m8/R/seed157,m8/R/seed905,m8/R/seed10157,m8/R/seed383,m8/R/seed10046,m8/R/seed981"


migrations=(mS M S R)
models=(sym asym oneRate threeRates)

# setup directories for xml files
for migration in ${migrations[@]}; do
for model in ${models[@]}; do
    dir="$output_dir/$migration/$model"
    mkdir -p $dir
    cp $input_xml $dir
done
done

for migration in ${migrations[@]}; do
for model in ${models[@]}; do
    if [[ "$model" == "sym" ]]; then
        symmetric="symmetric=true"
        spec="traitModelSpec=beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel"
        num_rates="numRates=55"
    elif [[ "$model" == "asym" ]]; then
        symmetric="symmetric=false"
        spec="traitModelSpec=beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel"
        num_rates="numRates=110"
    elif [[ "$model" == "oneRate" ]]; then
        symmetric="symmetric=true"
        spec="traitModelSpec=metastabayes.substitutionmodel.OneRateAllTissues"
        num_rates="numRates=1"
    elif [[ "$model" == "threeRates" ]]; then
        symmetric="symmetric=false"
        spec="traitModelSpec=metastabayes.substitutionmodel.ThreeRatesForSeedingRoutes"
        num_rates="numRates=3"
    fi

    echo -e "java -jar $metastabayes_jar -overwrite -working \
    -D "${input_names[$migration]}" \
    -D "$file_dir" \
    -D "$spec" \
    -D "$symmetric" \
    -D "$num_rates" \
    $output_dir/$migration/$model/proper_nested_sampling_machina_sims.xml > $output_dir/$migration/$model/${model}_terminal.log" >> $output_dir/parallel.txt
done
done
