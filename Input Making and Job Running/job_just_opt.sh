#!/bin/bash

#This script is meant to run an optimization job and a sinle-point unoptimized job given an .xyz file

#the .xyz file, without the ".xyz"
xyz=$1
#template qchem opt .in
opt=$2
#charge
charge=$3

#this version is written for Ulysses but only this line needs to be changed for telemachus (path to code)
python3 /work/ezraa/code/job_just_opt.py $xyz $opt $charge

sqthis -J "opt_${xyz}" -c 8 -t unlimited qchem.latest -nt 8 "opt_${xyz}.in" "opt_${xyz}.out"
