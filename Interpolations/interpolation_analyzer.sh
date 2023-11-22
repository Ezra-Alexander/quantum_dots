#!/bin/bash

#script to be run on laptop to collect data from a series of interpolated single points
#Inputs: Number of valence MO energies to extract, name of excel spreadsheet to write
#Should be run from the directory that contains all the frame subdirectories, and also contains nbas.txt

excel=$1 #name of excel file to write
nval=$2 #number of  MOs (starting from band edge) to plot #in the case of P4c, should be 3 per P
occ=$3 #either o for occupied orbitals or u for unoccupied

i=0
for d in */
do
	cd $d
	
	if [ $i != 0 ]
	then
		scp ../$excel .
	fi

	python3 ~/wormk/code/energies_2_excel.py $excel $nval $i $d $occ

	scp $excel ../

	cd ../
	i=$(( i + 1))
done