#!/bin/bash

#should take a qchem output file and run an aims dipole job on that geometry
#includes making geometry.in, editing control.in submit_aims.sh and submitting them

#To run this you need to have already copied a submit_aims.sh and qchem output to the directory, and need a qchem_plot.out and a control.in
#This will assume the same charge and multiplicity (and functional, basis etc.) as the control.in reference

#Goal - Take a qchem plot output file (already in the appropriate Telemachus directory) and generate and submit the FHIAims dipole matrix calculation
#Will require a bash script and a python script (or two)
#Needs to make the geometry file and edit the control and submit files

#the qchem output file
qchem=$1
#reference control.in
ref=$2

#run python script that makes geometry.in from qchem.out 
python3 ~/code/output_geom_qchem_2_aims.py $1

#rub python script that edits control.in and submit_aims.sh

python3 ~/code/control_aims.py $2 $1

#run the submit script
#sbatch submit_aims.sh