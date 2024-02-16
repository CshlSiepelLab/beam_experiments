#!/bin/bash

treefile=$1

figtree="/home/staklins/bin/FigTree_v1.4.4/lib/figtree.jar"



$figtree -graphic PDF $treefile
