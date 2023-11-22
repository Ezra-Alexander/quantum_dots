#!/bin/bash

#A script that runs a single point for each image of an interpolation, as performed by geodesic_interpolate 

images=$1 #the .xyz that contains all the images
ref=$2 #reference qchem single-point input, from which $rem is taken
charge=$3 #the charge of the system. Should be the same for all images
n_im=$4 #the number of images in the .xyz. geodesic_interpolate defaults to 17 (for some reason)

python3 ~/code/xyz_slicer.py $images #outputs xyz's 0.xyz to n_im.xyz

i=0
while  [ "$i" -lt "$n_im" ]; do
	mkdir "${i}th_frame"
	cd "${i}th_frame"

	scp "../${ref}" .
	mv "../${i}.xyz" .

	sh ~/code/job-forerunner.sh $i $ref $charge 

	cd ..
	i=$(( i + 1))
done
