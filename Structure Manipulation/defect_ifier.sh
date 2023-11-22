#!/bin/bash

#the goal here is to write a script capable of taking a .xyz file and a number of sets of atom indeces and create a subdirectory that includes subdirectories for each set, in which that set of atoms has been removed from the .xyz

#In other words, the script makes an arbitrary number of specified atom-removal defects, each packages neatly into different directories

#should not be run from a "defects" directory, as it makes said directory already
#only works on laptop because bash - run here then copy over
#which is to say it doesn't run them

#the original .xyz. Assumed to begin with "optzd_"
xyz=$1

shift 1

#the following arguments should be structured as follows: n1 i1 i2 ... in1 n2 ... nm  where the nj are integers specifying the number of atoms to be removed in the jth defect being formed, m is the total number of defects and i1 ... inj specify the overall indeces of each of the atoms removed in the jth defect

#make the defects subdirectory. It is made this way so that it is easy to copy the entire thing to the cluster
mkdir defects
cd defects

#now we begin the loop that 
while [[ $# -gt 0 ]]; do

	num=$1

	#echo "${@:2:($num)}"

	#Now what should I name the subdirectories for each defect ...
	name=""
	for i in ${@:2:($num)}; do
		name="${name}${i}_"
	done
	name=${name%_}

	mkdir $name
	cd $name

	scp ../../$xyz .

	python3 ~/wormk/code/defect_ifier.py $xyz ${@:2:($num)}

	shift_val=$(( 1+$num ))

	shift $shift_val

	cd ../

	if [[ $# -gt 0 ]]; then
		if [[ $1 -gt $# ]]; then
			echo "You didn't format your input correctly, dingus"
			break
		fi
	fi

done

cd ../
