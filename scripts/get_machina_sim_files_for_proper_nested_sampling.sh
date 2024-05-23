#!/bin/bash

### Basic bash scripting to find and format all machina sim files to run nested sampling across a dataset

# user input desired working dir
working_dir="machina_sims_proper_nested_sampling_5_23_24"
mkdir -p $working_dir

# find all sim dirs based on existing location
sim_dirs=$(find machina_data/sims/*/ -type d -name seed*)

# set empty string for beast command line string
sim_names_string=""

i=0

# main file processing steps here
for dir in $sim_dirs; do
# make new name id for sim to condense dir string
m5_or_m8=$(echo $dir | cut -d'/' -f3 | cut -d'_' -f2)
seed=$(echo $dir | cut -d'/' -f4)
new_name=$m5_or_m8"_"$seed

# get newick and tissues files for each sim
newick=$dir/*_unlabeled_true_tree.nwk
tissues=$dir/*_tissues.tsv

# copy newick and tissues files to working dir with new name
new_newick=$working_dir/$new_name".nwk"
new_tissues=$working_dir/$new_name".tsv"
cp $newick $new_newick
cp $tissues $new_tissues

# call feast file io reformatting for each sim
./scripts/setup_feast_io_files_nested_sampling.sh $new_newick $new_tissues

# add to string input for beast command line of all sim names
if [[ "$sim_names_string" == "" ]]; then
    sim_names_string+="$new_name"
else
    sim_names_string+=",$new_name"
fi

# count the number of sims
i=$(( i+1 ))
done

echo "inputNames=\"$sim_names_string\"" > $working_dir/commandline_input.txt
echo "$i sims processed"
