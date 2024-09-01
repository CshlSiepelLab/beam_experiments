#!/bin/bash

### This script takes in a single xml and determines the subchainlength for nested sampling in beast with one active particle as recommended by the developers
### Could probably implement a more efficient version of this script, but need to consider when smaller chain lengths satisfy the diff threshold and then larger chains do not do to stochastic nature (could maybe aggregarte across multithreaded runs for multiple nested sampling results for each subChainLength)

###########
# Info from BEAST nested sampling FAQs on subChainLength:
# "NS works in theory if and only if the points generated at each iteration are independent. 
# If you already did an MCMC run and know the effective sample size (ESS) for each parameter, 
# to be sure every parameter in every sample is independent you can take the length of the MCMC 
# run divided by the smallest ESS as sub-chain length. This tend to result in quite large sub-chain 
# lengths.
# In practice, we can get away much smaller sub-chain lengths, which you can verify by running multiple 
# NS analysis with increasing sub-chain lengths. If the ML and SD estimates do not substantially differ, 
# (you want the difference between ML1 and ML2 to be at most 2*sqrt(SD1*SD1+SD2*SD2)) you know the shorter 
# sub-chain length was sufficient."
###########


if [[ $# -eq 0 ]] ; then
    echo "Usage: determine_subchainlength_nested_sampling.sh --xml <xml filepath (str)>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -x|--xml) xml="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: determine_subchainlength_nested_sampling.sh --xml <xml filepath (str)>"; exit 1 ;;
    esac
    shift
done

xml_filename=$(basename $xml)
beast_path=$(which beast)

# set subchain lengths to test
test=(100 500 1000 2500 5000 10000 20000)

# set empty array to obtain ML and SD results
ml_values=()
sd_values=()

for subchainlength in ${test[@]};
do

echo "Running subChainLength of $subchainlength"

# make temporary xml copy to prep for nested sampling
mkdir temp_dir
temp_xml="temp_dir/temp_${xml_filename}"
cp $xml $temp_xml

# replace necessary portions of the xml to run nested sampling
existing_mcmc="<run id=\"mcmc\" spec=\"MCMC\""
active_particles=1
# can setup nested sampling with or without autoSubChainLength which may be unstable
replace_mcmc="<run id=\"mcmc\" spec=\"beast.gss.NS\" chainLength=\"1000000\" particleCount=\"$active_particles\" subChainLength=\"$subchainlength\" epsilon=\"1e-13\">"
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$temp_xml"

existing_logger='<logger id="tracelog" spec="Logger'
replace_logger='<logger id="tracelog" spec="NSLogger'
sed -i "s|$existing_logger|$replace_logger|" "$temp_xml"

# run beast for temp xml
output_xml="temp_dir/xml1_marginal_likelihood_run.txt"
$beast_path -overwrite -working "$temp_xml" > "$output_xml"

# obtain ML and SD from the run output
ml_line_xml=$(grep "Marginal likelihood:" "$output_xml" | grep -vE "subsample|bootstrap")
ml=$(echo "$ml_line_xml" | cut -d ' ' -f 3 | cut -d '(' -f 1)
sd=$(echo "$ml_line_xml" | cut -d ' ' -f 3 | cut -d '(' -f 2 | cut -d ')' -f 1)

# add ML and SD to array to track results
ml_values+=("$ml")
sd_values+=("$sd")

# remove temp files after each run
rm -r temp_dir

done


# determine if differences between subChainLength ML estimates exceed 2*sqrt(sd1^(2)+sd2^(2)) to find at what subChainLength the threshold is not exceeded anymore
diffs=()

for ((i=1; i<${#ml_values[@]}; i++)); do

prev_index=$(( $i-1 ))
ml1="${ml_values[prev_index]}"
sd1="${sd_values[prev_index]}"
ml2="${ml_values[i]}"
sd2="${sd_values[i]}"

log_bf=$(echo "$ml1 - $ml2" | bc -l)
abs_log_bf=$(echo "if ($log_bf < 0) -($log_bf) else $log_bf" | bc -l)
diff_threshold=$(echo "2 * sqrt(($sd1^2) + ($sd2^2))" | bc -l)

echo "subChainLength: ${test[prev_index]} to ${test[i]}     Absolute value of log(BF): ${abs_log_bf}      Threshold: ${diff_threshold}"

if (( $(echo "$abs_log_bf < $diff_threshold" | bc -l) )); then
    diffs+=("1")
else
    diffs+=("0")
fi

done

# find the last ocurence of log(BF) less than threshold where there are no false comparisons to the right of it, indicating that stability has been reached
barrier_index=-1

for ((i=0; i<${#diffs[@]}; i++)); do
    # Check if the current element is 1
    if [[ "${diffs[i]}" -eq 1 ]]; then
        # Check if there are no more 0 values to the right
        if [[ ! "${diffs[@]:$i}" =~ 0 ]]; then
            # Set the barrier_index and break out of the loop
            barrier_index=$i
            break
        fi
    fi
done

echo ${diffs[@]}
echo $barrier_index

if [ "$barrier_index" -eq -1 ]; then
    echo "Error: ML stability not reached, so please test larger subChainLengths!"
else
    step_size="${test[barrier_index]}"
fi

echo "Use a step size of ${step_size}"
