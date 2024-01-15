#!/bin/bash

#This script is meant to run an optimization job and a sinle-point unoptimized job given an .xyz file

#the .xyz file. Note that the file name is assumed to be unopt_DESCRIPTOR.xyz
xyz=$1
#template qchem opt .in
opt=$2
#charge
charge=$3
#number of corees
cores=$4
#priority. If on Ulysses, leave blank
priority=$5

job_name=${xyz%.*}
job_name=${job_name#*_}

python3 ~/code/job_just_opt.py $xyz $opt $charge

#note that on Ulysses the bash version requires if [ ]; instead
if [[ -z $priority ]]; then
	sqthis -J "opt_${job_name}" -c $cores -t unlimited qchem.latest -nt $cores "opt_${job_name}.in" "opt_${job_name}.out"
else
	sqthis -J "opt_${job_name}" -c $cores -t unlimited -p $priority qchem.latest -nt $cores "opt_${job_name}.in" "opt_${job_name}.out"
fi

