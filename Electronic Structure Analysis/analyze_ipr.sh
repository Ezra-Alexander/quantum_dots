#!/bin/bash

#this script will perform the on-cluster portions of ipr analysis - making an .xyz file, finding the rows and basis functions, and running get_low_perorb and ipr_orb
#for now it will not loop over multiple directories - just the one
#before running, I may need to copy the ipr output from ulysses to telemachus
#after running, I need to copy the ultimate output - low_orb_ipr.csv - to my laptop, in the respective pdos directory
#configured for telemachus

#the ipr output file
ipr=$1
#the plot input file
plot=$2

#make .xyz file
python3 ~/wormk/code/make_xyz.py $plot

#get rows and basis functions, saved to rows_basis.txt
#note that this doesn't work for rows if the structure has H in it
python3 ~/wormk/code/rows_and_basis_functions_er.py $ipr

#run the new version of get_low_perorb that accepts a txt file instead of integers
python3 ~/wormk/code/get_low_perorb_ezra.py $ipr optzd* rows_basis.txt

#run ipr_orb.py
python3 ~/wormk/code/ipr_orb.py optzd* low_orb.npy
