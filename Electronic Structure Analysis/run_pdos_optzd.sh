#!/bin/bash

#run from outside directory of job(s). First argument is path on ulysses to directory of job(s), subsequent arguments are actual job directories of interest

i=0

for dir in "$@"
do
	if [ "$i" -eq 0 ]
	then
		i=$(($i+1))
		continue
	fi

	cd $dir
	#currently configured for optzd single points
	mkdir "pdos_optzd"
	cd "pdos_optzd"
	
	#needs to be changed if using telemachus
	cluster=ezraa@ulysses:
	#assumes name of directory on cluster matches name of directory here
	path="$cluster$1"
	path_53="${path}/${dir}/scratch*optzd*/53.npz"
	path_nbas="${path}/${dir}/scratch*optzd*/nbas.txt"
	path_in="${path}/${dir}/*optzd*.in"


	scp $path_53 .
	scp $path_nbas .
	scp $path_in .

	#we assume there is no .xyz file for the optimized structure yet and make it ourselves

	python3 ~/wormk/code/make_xyz.py *optzd*.in

	echo $dir 

	#note that this works only for the pdos code for InP dots
	python3 /Users/ezraalexander/wormk/code_from_desktop/code/pdos_inp.py *.xyz nbas.txt 53.npz


	cd ../..

done
