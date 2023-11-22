#!/bin/bash

#Developed with F to hydroxide substitutions in mind but it would be nice to make this more general
#the idea is to take a dot, choose a type of substitution defect (e.g. F- to OH-)
#then create some specified number of different random versions of this defect, then run opts for each

#relies on qchem helper
#supports hydroxides and oxides
#otherwise assumes simple atomic substitution with no bond length adjustment

#file to get the optimized structure from. Takes .xyz directly, opt.out, or plot_optzd.in
xyz=$1

#reference opt.in
ref=$2

#atom to be replaced
orig=$3

#specied to replace with. For species with more than one atom, specify both in one word. i.e. "OH" for hydroxide
replace=$4

#number of random substitutions to make
nsub=$5

#charge of your dot after substitution (there's a way to write this code so I don't need to specify this, but this is easier)
charge=$6

#the meat of the script. calls python script to make the 5 random .xyzs
#path needs to be changed for running on different clusters
python3 ~/code/random_substitutioner.py $xyz $orig $replace $nsub

#now make the subdirectories and run the jobs
i=1
while  [ "$i" -le "$nsub" ]; do
	mkdir "${i}th_sub"
	cd "${i}th_sub"

	scp "../${ref}" .
	mv "../${i}.xyz" .

	#again, be careful with your path
	sh ~/code/job_just_opt.sh $i $ref $charge 

	cd ..
	i=$(( i + 1))
done