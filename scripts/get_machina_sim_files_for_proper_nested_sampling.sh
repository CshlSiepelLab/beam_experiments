#!/bin/bash

### Basic bash scripting to find and format all machina sim files to run nested sampling across a dataset

# user input desired working dir
working_dir="machina_sims_proper_nested_sampling_5_29_24"
mkdir -p $working_dir

# copy machina sims to working dir
cp -r machina_data/sims/* $working_dir/

# remove unnecessary files
find $working_dir/ -type f -name "*_labeled_true_tree.nwk" -delete
find $working_dir/ -type f -name "*.tree" -delete
find $working_dir/ -type f -name "*.vertex.labeling" -delete

# rename files to seed only .nwk and .tsv
files=$(find $working_dir/ -type f -name "*_unlabeled_true_tree.nwk")
for file in $files; do dir=$(dirname $file); echo $dir; seed=$(echo $file | cut -d'_' -f 9); echo $seed; mv $file $dir/$seed.nwk; done
files=$(find $working_dir/ -type f -name "*.tsv")
for file in $files; do dir=$(dirname $file); echo $dir; seed=$(echo $file | cut -d'_' -f 9); echo $seed; mv $file $dir/$seed.tsv; done

# call feast file io reformatting for each sim
files=$(find $working_dir/ -type f -name "*.nwk")
for file in $files; do
nwk=$file
tsv=${file//.nwk/.tsv}
./scripts/setup_feast_io_files_nested_sampling.sh $nwk $tsv
done

# fix bug in tips csv
files=$(find $working_dir -type f -name *_tips_only.csv)
for file in $files; do
sed -i 's/ /\n/g' $file
done

# add to string input for beast command line of all sim names
files=$(find $working_dir/ -type f -name "*_reformatted.nwk")
sim_string=""
for file in $files; do
sim_name=$(echo $file | cut -d'_' -f 8 | cut -d'/' -f 2-4)
sim_string+="$sim_name,"
done

