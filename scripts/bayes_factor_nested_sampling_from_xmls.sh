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
existing_mcmc='<run id="mcmc" spec="MCMC"'
replace_mcmc='<run id="mcmc" spec="beast.gss.NS" chainLength="250000" particleCount="1" subChainLength="20000" epsilon="1e-13" autoSubChainLength="true" paramCountFactor="10">'
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$cp_xml1"
sed -i "s|$existing_mcmc.*|$replace_mcmc|" "$cp_xml2"

existing_logger='<logger id="tracelog" spec="Logger'
replace_logger='<logger id="tracelog" spec="NSLogger'
sed -i "s|$existing_logger|$replace_logger|" "$cp_xml1"
sed -i "s|$existing_logger|$replace_logger|" "$cp_xml2"

# run nested sampling in BEAST for both xml files
beast_path=$(which beast)
$beast_path -overwrite -working $cp_xml1 > ${xml1_dir}/xml1_marginal_likelihood_run.txt
$beast_path -overwrite -working $cp_xml2 > ${xml2_dir}/xml2_marginal_likelihood_run.txt




