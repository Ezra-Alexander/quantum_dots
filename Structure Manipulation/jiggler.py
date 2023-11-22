import numpy as np
import sys
import math
from qchem_helper import *
import random

#a simple script that takes a qchem opt input file and "jiggles" the structure by moving all atoms a random (small) amount

#originally for testing GaP defect convergence

opt_in = sys.argv[1] #the qchem opt.in you're changing
opt_out=sys.argv[2] #the name of the opt.in you want to write

wiggle=0.1 #manually set parameter. Twice the max distance an atom is allowed to move, in Angstroms

atoms,coords = get_geom_io(opt_in)
rem,sp=my_get_rem_sp(opt_in)

new_rem=""
for i,line in enumerate(rem):
	new_rem=new_rem+line

for i,atom in enumerate(coords):
	for j,coord in enumerate(atom):
		perturbation=random.randrange(0,100)
		perturbation=(perturbation/100)-0.5
		perturbation=perturbation*wiggle
		atom[j]=coord+perturbation
		

write_input(opt_out,atoms,coords,new_rem,sp)

