#!/bin/bash

#This script originally ran an optimization job and a sinle-point unoptimized job given an .xyz file
#I'm going to modify it to work with just the single point

#the .xyz file, without the ".xyz"
xyz=$1
#template qchem opt .in
#opt=$2
#template qchem sp .in
sp=$2
#charge
charge=$3

#this version is written for telemachus but only this line needs to be changed for ulysses (path to code)
python3 ~/code/job_forerunner.py $xyz $sp $charge $subscrip

#sqthis -J "opt_${xyz}" -c 8 -t unlimited qchem.latest -nt 8 "opt_${xyz}.in" "opt_${xyz}.out"

sbatch submit_plot.sh
