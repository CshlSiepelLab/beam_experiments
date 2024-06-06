#!/bin/bash

nonM1uniformProb=0.05263
nonM1lastProb=0.05266
uniformProb=0.02778
lastProb=0.02774
m1prob=0.5
zeroProb=0.0
dimension=20

for ((i=1; i<=dimension; i++))
do
    for ((j=1; j<=dimension; j++))
    do
        # if last row first element
        if [ $i -eq $dimension ] && [ $j -eq 1 ]; then
            printf "$lastProb"
        # # if second row diagonal element
        # elif [ $i -eq 2 ] && [ $i -eq $j ]; then
        #     printf "1.0"
        # if diagonal element
        elif [ $i -eq $j ]; then
            printf "$zeroProb"
        # # if second row
        # elif [ $i -eq 2 ]; then
        #     printf "$zeroProb"
        # if second row last element
        elif [ $i -eq 2 ] && [ $j -eq $dimension ]; then
            printf "$nonM1lastProb"
        elif [ $i -eq 2 ]; then
            printf "$nonM1uniformProb"
        # if second row
        elif [ $i -eq 2 ]; then
            printf "$zeroProb"
        # if second element in any row
        elif [ $j -eq 2 ]; then
            printf "$m1prob"
        # if last element in any row
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