#!/bin/bash
#useage: script is run with ./analysis.sh -b basepath -i inputfilepath
#where the inputfile is the raw out from CP2K
OPTIND=1

input_path=""
base_path=""

while getopts "b:i:" opt; do
    case "$opt" in
    b)  base_path=$OPTARG
        ;;
    i)  input_path=$OPTARG
        ;;
    esac
done

grep "forces " $base_path$input_path  | cut -d ' ' -f2- | sort > $base_path$input_path".temp"

python condenseforces.py $base_path$input_path".temp" $base_path$input_path".final"

rm $base_path$input_path".temp"
