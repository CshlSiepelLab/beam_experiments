#!/bin/bash

m5_seeding_route_tsv="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/machina_data/machina_sims_m5_seeding_routes_key.tsv"
m8_seeding_route_tsv="/grid/siepel/home_norepl/staklins/bayesian_phylogenetic_metastasis/machina_data/machina_sims_m8_seeding_routes_key.tsv"

# Read m5_seeding_route_tsv into a dictionary
unset m5_seeding_route_dict
declare -A m5_seeding_route_dict
while IFS=$'\t' read -r key value; do
    m5_seeding_route_dict["$key"]+=" m5_$value"
done < "$m5_seeding_route_tsv"

# Read m8_seeding_route_tsv into a dictionary
unset m8_seeding_route_dict
declare -A m8_seeding_route_dict
while IFS=$'\t' read -r key value; do
    m8_seeding_route_dict["$key"]+=" m8_$value"
done < "$m8_seeding_route_tsv"

# make combined dictionary
unset combined_dict
declare -A combined_dict
for key in "${!m5_seeding_route_dict[@]}"; do
    combined_dict["$key"]="${m5_seeding_route_dict[$key]}${m8_seeding_route_dict[$key]}"
done

# sort files
for key in "${!combined_dict[@]}"; do
    for seed in ${combined_dict[$key]}; do
        mv machina_sims_proper_nested_sampling_5_29_24/$seed[._]* machina_sims_proper_nested_sampling_5_29_24/$key/
    done
done

# print comma sep string for each key
for key in "${!combined_dict[@]}"; do
    if [[ $key == "Subdirectory" ]]; then
        continue
    fi
    echo $key
    echo ${combined_dict[$key]} | tr ' ' ','
    # echo ${combined_dict[$key]} | tr ' ' ',' | sed "s/,/,$key\//g" | awk -v key="$key/" '{print key$0}'
done
