#!/bin/bash

#for scratch directories with "unopt" in their name
#nbas.txt needs to be in same directory

for dir in "$@"
do
	cd $dir
	cd scratch*unopt*
	scp ../../nbas.txt .
	sh /work/ezraa/code/scratch_to_npy.sh
	echo Done with $dir
	echo
	cd ../../
done
