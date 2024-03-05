#!/bin/bash

# A script that takes a .xyz
# And runs a q-chem single point that includes some combination of:
#	saved scratch for PDOS
#   plots of specified orbitals
#   unoccupied lowdin populations for the IPR
# The focus is to make this as easy to run as possible
# Also automatically 0-centers your coordinates just for funsies

# Assumptions:  State in question is a singlet
#				That the joint plot/ipr job's scratch will save the way I want it to (it does)
#				That we will be using PBE0/def2-SVP
#				Everything is in its most common oxidation state (+1 for H)
#				.xyz file name is of the form optzd_DESCRIPTOR.xyz
#				That we want to use rca_diis as our scf algorithm
#				That our dot is small enough to fit in a (-20,20) grid

# as of now, only veryhigh priority has been implemented but the framework is there

xyz=$1 #the .xyz file. Hopefully the only file input
cores=$2 #the number of cores you want the calculation to run on
plot_min=$3 #the lowest energy MO you want to plot. 0 for HOMO, 1 for HOMO-1, etc
plot_max=$4 #the highest energy MO you want to plot. 0 for LUMO, 1 for LUMO+1, etc
mode=$5 #just a plot ("plot"), just an ipr ("ipr"), or the ipr+plot ("both")
mem_static=$6 #the mem_static value in qchem's $rem, in MB
priority=$7 #the priority for the job. If on Ulysses, use 'ulysses'. priotity automatically sets mem_total

python3 ~/code/xyz_to_single_point.py $xyz $cores $plot_min $plot_max $mode $mem_static $priority

sbatch submit_plot.sh