#!/bin/bash

### This script uses nested sampling in BEAST to compute marginal likelihood given a model with data in an xml file

# set defaults for if user input does not specify these conditions
# start with low number of active particles such as 1 and then increase only if SD of each is not less than 2 so the diff threshold is too high. Can compute the necessary number of particles with N=H/(SD*SD) where H is the information content and SD the desired standard deviation, reccomended to be 2.
active_particles=1
# can use another script to test different subChainLengths to choose that which is the minimum requried to ensure ML and SD stability which indicates indepedent sampling of points. Values in the range of 5000-20000 seem to be the most common. Can also empirically determine based on seperate normal MCMC run where subChainLength = len MCMC / smallest ESS of all parameters, but this will give an unneccessarily high value.
sub_chain_length=10000

if [[ $# -eq 0 ]] ; then
    echo "Usage: nested_sampling_marginal_likelihood_from_xml.sh --xml <xml filepath (str)> --dir <working directory path (str)> --active_particles <(int)> --sub_chain_len <(int)>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -x|--xml) xml1="$2"; shift ;;
        -d|--dir) dir="$2"; shift ;;
        -p|--active_particles) active_particles="$2"; shift ;;
        -l|--sub_chain_length) sub_chain_length="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: nested_sampling_marginal_likelihood_from_xml.sh --xml <xml filepath (str)> --dir <working directory path (str)>"; exit 1 ;;
    esac
    shift
done


# make independent working directories for the nested sampling runs to avoid interference
xml1_dir="${dir}/xml1"

mkdir -p $xml1_dir

# copy the provided xml files to the workign directories to make modifications for a nested sampling run
cp $xml1 $xml1_dir
xml1_file=$(basename "$xml1")
cp_xml1="${xml1_dir}/${xml1_file}"

# replace necessary portions of the xml to run nested sampling
existing_mcmc="<run id=\"mcmc\" spec=\"MCMC\""
replace_mcmc="<run id=\"mcmc\" spec=\"beast.gss.NS\" chainLength=\"1000000\" particleCount=\"$active_particles\" subChainLength=\"$sub_chain_length\" epsilon=\"1e-13\">"
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$cp_xml1"

existing_logger='<logger id="tracelog" spec="Logger'
replace_logger='<logger id="tracelog" spec="NSLogger'
sed -i "s|$existing_logger|$replace_logger|" "$cp_xml1"

# run nested sampling in BEAST for both xml files
beast_path=$(which beast)
output_xml1="${xml1_dir}/xml1_marginal_likelihood_run.txt"

# run xml nested sampling
$beast_path -overwrite -working $cp_xml1 > $output_xml1

