#!/bin/bash

uniformProb=0.02778
lastProb=0.02774
m2prob=0.5
zeroProb=0.0
dimension=20

for ((i=1; i<=dimension; i++))
do
    for ((j=1; j<=dimension; j++))
    do
        if [ $i -eq $j ]; then
            printf "$zeroProb"
        elif [ $j -eq 2 ]; then
            printf "$m2prob"
        elif [ $j -eq $dimension ]; then
            printf "$lastProb"
        else
            printf "$uniformProb"
        fi

        if [ $j -ne $dimension ]; then
            printf ","
        fi
    done
    printf "\n"
done