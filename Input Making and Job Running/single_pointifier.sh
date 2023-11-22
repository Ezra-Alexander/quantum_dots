#!/bin/bash

#this bash script runs a single point for a converged opt
#should call the titular python script, which makes the single point input file and the submit script
#it then runs the submit script
#written for ulysses

#input is the name of the input file for the opt you already ran, not the name of the sp input you want to make
input=$1
output=$2

python /work/ezraa/code/single_pointifier.py $input $output

sbatch new_submit.sh