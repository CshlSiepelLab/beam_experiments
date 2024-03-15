#!/bin/bash

### THis script takes two output files with nested sampling marginal likelihood from BEAST reported and calculates bayes factor between them for interpretation

if [[ $# -eq 0 ]] ; then
    echo "Usage: bayes_factor_from_marginals.sh --ml1 <marginal likelihood output 1 filepath (str)> --ml2 <marginal likelihood output 2 filepath (str)> --dir <working directory path (str)>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -x|--ml1) output_xml1="$2"; shift ;;
        -y|--ml2) output_xml2="$2"; shift ;;
        -d|--dir) dir="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: bayes_factor_from_marginals.sh --ml1 <marginal likelihood output 1 filepath (str)> --ml2 <marginal likelihood output 2 filepath (str)> --dir <working directory path (str)>"; exit 1 ;;
    esac
    shift
done



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

if (( $(echo "$abs_log_bf <= 1.1" | bc -l) )); then
    intepretation="${model} is preferred, but the difference is hardly worth mentioning."
    echo $interpretation
elif (( $(echo "$abs_log_bf <= 3" | bc -l) )); then
    interpretation="${model} is preferred with positive support."
    echo $interpretation
elif (( $(echo "$abs_log_bf <= 5" | bc -l) )); then
    interpretation="${model} is preferred with strong support."
    echo $interpretation
elif (( $(echo "$abs_log_bf > 5" | bc -l) )); then
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
