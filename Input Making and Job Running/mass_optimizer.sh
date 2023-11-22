#!/bin/sh

#should take a number of .xyz files and run job_just_opt.sh for all of them
#assumes charge 0

ref=$1

for .xyz in "$@"
do

	if [ "$i" -eq 0 ]
	then
		i=$(($i+1))
		continue
	fi

	#configured for Ulysses
	sh /work/ezraa/code/job_just_opt.sh $xyz $ref 0

done