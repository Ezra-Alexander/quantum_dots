#!/bin/bash

#first command line argument is the path on ulysses to the directory holding the individual directories for each job, and the following arguments are each of those jobs
#currently configured for ulysses

i=0

for dir in "$@"
do
	if [ "$i" -eq 0 ]
	then
		i=$(($i+1))
		continue
	fi

	cd $dir
	
	#needs to be changed if using ulysses
	cluster=ezraa@ulysses:
	#assumes name of directory on cluster matches name of directory here
	path="$cluster$1"
	#currently configured for unopt single points
	final_path="${path}/${dir}/*.plots"


	scp -r $final_path .

	cd ..

done
