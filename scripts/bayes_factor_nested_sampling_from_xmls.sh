#!/bin/bash

### This script takes in two xml files and modifies them as necessary to setup xml files for nested sampling analysis to get marginal likelihoods, compute a bayes factor, and select between the models

if [[ $# -eq 0 ]] ; then
    echo "Usage: bayes_factor_nested_sampling_from_xmls.sh --xml1 <xml filepath (str)> --xml2 <xml filepath (str)> --dir <working directory path (str)>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -x|--xml1) xml1="$2"; shift ;;
        -y|--xml2) xml2="$2"; shift ;;
        -d|--dir) dir="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: Usage: bayes_factor_nested_sampling_from_xmls.sh --xml1 <xml filepath (str)> --xml2 <xml filepath (str)> --dir <working directory path (str)>"; exit 1 ;;
    esac
    shift
done


# make independent working directories for the nested sampling runs to avoid interference
xml1_dir="${dir}/xml1"
xml2_dir="${dir}/xml2"

mkdir -p $xml1_dir
mkdir -p $xml2_dir

# copy the provided xml files to the workign directories to make modifications for a nested sampling run
cp $xml1 $xml1_dir
cp $xml2 $xml2_dir
xml1_file=$(basename "$xml1")
xml2_file=$(basename "$xml2")
cp_xml1="${xml1_dir}/${xml1_file}"
cp_xml2="${xml2_dir}/${xml2_file}"

# replace necessary portions of the xml to run nested sampling
existing_mcmc="<run id=\"mcmc\" spec=\"MCMC\""
# start with low number of active particles such as 1 and then increase only if SD of each is not less than 2 so the diff threshold is too high. Can compute the necessary number of particles with N=H/(SD*SD) where H is the information content and SD the desired standard deviation, reccomended to be 2.
active_particles=1
# can use another script to test different subChainLengths to choose that which is the minimum requried to ensure ML and SD stability which indicates indepedent sampling of points. Values in the range of 5000-20000 seem to be the most common. Can also empirically determine based on seperate normal MCMC run where subChainLength = len MCMC / smallest ESS of all parameters, but this will give an unneccessarily high value.
subchainlength=10000
replace_mcmc="<run id=\"mcmc\" spec=\"beast.gss.NS\" chainLength=\"1000000\" particleCount=\"$active_particles\" subChainLength=\"$subchainlength\" epsilon=\"1e-13\">"
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$cp_xml1"
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$cp_xml2"

existing_logger='<logger id="tracelog" spec="Logger'
replace_logger='<logger id="tracelog" spec="NSLogger'
sed -i "s|$existing_logger|$replace_logger|" "$cp_xml1"
sed -i "s|$existing_logger|$replace_logger|" "$cp_xml2"

# run nested sampling in BEAST for both xml files
beast_path=$(which beast)
output_xml1="${xml1_dir}/xml1_marginal_likelihood_run.txt"
output_xml2="${xml2_dir}/xml2_marginal_likelihood_run.txt"

# run both xml files in parallel
cmd1="$beast_path -overwrite -working '$cp_xml1' > '$output_xml1'"
cmd2="$beast_path -overwrite -working '$cp_xml2' > '$output_xml2'"
commands=("$cmd1" "$cmd2")

# setup file with 2 individual commands to run in parallel
for command in "${commands[@]}"
do
  echo "${command}" >> "${dir}/parallel.txt"
done

# run 2 individual commands in parallel
parallel -j 2 < "${dir}/parallel.txt"
rm "${dir}/parallel.txt"

# grab marginal likleihood and SD estimates from the plain output
ml_line_xml1=$(grep "Marginal likelihood:" "$output_xml1" | grep -vE "subsample|bootstrap")
ml_line_xml2=$(grep "Marginal likelihood:" "$output_xml2" | grep -vE "subsample|bootstrap")

# process output line to get ML and SD for each xml run
ml1=$(echo "$ml_line_xml1" | cut -d ' ' -f 3 | cut -d '(' -f 1)
sd1=$(echo "$ml_line_xml1" | cut -d ' ' -f 3 | cut -d '(' -f 2 | cut -d ')' -f 1)
ml2=$(echo "$ml_line_xml2" | cut -d ' ' -f 3 | cut -d '(' -f 1)
sd2=$(echo "$ml_line_xml2" | cut -d ' ' -f 3 | cut -d '(' -f 2 | cut -d ')' -f 1)

# compute Bayes factor and difference threshold from ML and SD estimates for each XML
log_bf=$(echo "$ml1 - $ml2" | bc -l)
abs_log_bf=$(echo "if ($log_bf < 0) -($log_bf) else $log_bf" | bc -l)
diff_threshold=$(echo "2 * sqrt(($sd1^2) + ($sd2^2))" | bc -l)

# determine if the number of particles is enough
if (( $(echo "$abs_log_bf < $diff_threshold" | bc -l) )); then
    echo "Bayes factor comparison failed because the Log(BF) of $abs_log_bf was not greater than $diff_threshold. Need to repeat with more active particles then ${active_particles} to reduce the difference threshold to detect fine differences, or otherwise conclude that the models cannot be distinguished."
    exit
fi

# decide if model 1 or model 2 is favored
if (( $(echo "$log_bf > 0" | bc -l) )); then
    model="Model1 (xml1)"
elif (( $(echo "$log_bf < 0" | bc -l) )); then
    model="Model1 (xml1)"
fi

# interpret log Bayes factor where a positive value support model1 and a negative value supports model2
# echo "log(BF) = log(ML1) - log(ML2), where BEAST nested sampling automatically reports log(ML) values"
echo "Log(BF) is ${log_bf}"

if (( $(echo "$abs_log_bf <= 0.5" | bc -l) )); then
    intepretation="${model} is preferred, but the difference is hardly worth mentioning."
    echo $interpretation
elif (( $(echo "$abs_log_bf <= 1.3" | bc -l) )); then
    interpretation="${model} is preferred with positive support."
    echo $interpretation
elif (( $(echo "$abs_log_bf <= 2.2" | bc -l) )); then
    interpretation="${model} is preferred with strong support."
    echo $interpretation
elif (( $(echo "$abs_log_bf > 2.2" | bc -l) )); then
    interpretation="${model} is preferred with overwhelming support."
    echo $interpretation
fi

# write results to an output file
outputfile="${dir}/bayes_factor_results.txt"

echo "Log(ML1) is ${ml1}" > $outputfile
echo "SD of ML1 is ${sd1}" >> $outputfile
echo "Log(ML2) is ${ml2}" >> $outputfile
echo "SD of ML2 is ${sd2}" >> $outputfile
echo "Log(BF) is ${log_bf}" >> $outputfile
echo $interpretation >> $outputfile
