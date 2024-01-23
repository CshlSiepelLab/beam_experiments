#!/bin/bash

# This script runs tidetree.jar from an input xml file
# Output will be in the working directory from which this script is called

if [[ $# -eq 0 ]] ; then
    echo "Usage: run_tidetree.sh --xml <tidetree xml filepath (str)>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -x|--xml) XML="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: run_tidetree.sh --xml <tidetree xml filepath (str)>"; exit 1 ;;
    esac
    shift
done

TIDETREE_JAR="../../tidetree/bin/tidetree.jar"

if [ ! -f $tidetree_jar ]
then
    echo "Script ../../tidetree/bin/tidetree.jar not found. Please install tidetree repo upstream of this repo or correct path in this file to tidetree.jar. Exiting!"
    exit
fi



java -jar $TIDETREE_JAR $XML


