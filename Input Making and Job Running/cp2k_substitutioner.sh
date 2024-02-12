#!/bin/bash

#a script for rapidly running a number of geometry optimizations in cp2k
#for a number of unique randomly substituted versions of an initital structure
#makes X copies each with Y unique random substitutions

#file to get the optimized structure from. Takes .xyz directly
xyz=$1

#a reference cp2k opt.in. Only the geometry is changed
ref_opt=$2

#a reference cp2k submit script. only file names are changed
ref_sub=$3

#element to replace
orig=$4

#element to add
new=$5

#number of structures to make
nstruct=$6

#number of substitutions per structure. If less than 1, treated as a fraction
nsub=$7

#python script to make all the random .xyzs in subdirectories
python3 ~/wormk/code/mass_substitutioner.py $xyz $orig $new $nstruct $nsub

#now, loop over all subdirectories and run the opts
for d in */
do
	cd $d
	
	python3 ~/wormk/code/cp2k_opt_maker.py unopt_*.xyz ../$ref_opt ../$ref_sub
	echo sbatch run_cp2k.sh

	cd ../

done