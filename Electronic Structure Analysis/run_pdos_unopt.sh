#!/bin/bash

i=0

for dir in "$@"
do
	if [ "$i" -eq 0 ]
	then
		i=$(($i+1))
		continue
	fi

	cd $dir
	#currently configured for unopt single points
	mkdir "pdos_unopt"
	cd "pdos_unopt"
	
	#needs to be changed if using ulysses
	cluster=ezraa@telemachus:
	#assumes name of directory on cluster matches name of directory here
	path="$cluster$1"
	path_53="${path}/${dir}/scratch*unopt*/53.npz"
	path_nbas="${path}/${dir}/scratch*unopt*/nbas.txt"


	scp $path_53 .
	scp $path_nbas .

	#for unopt there should be a .xyz file in the directory on my computer

	scp ../*.xyz .

	echo $dir 

	#note that this works only for the pdos code for InP dots
	python3 /Users/ezraalexander/wormk/code_from_desktop/code/pdos_inp.py *.xyz nbas.txt 53.npz


	cd ../..

done