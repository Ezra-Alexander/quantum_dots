import numpy as np
import sys
import math
from geom_helper import *

#quick script for figuring out problems with structures

xyz=sys.argv[1]

#read xyzs
coords, atoms = read_input_xyz(xyz)


def_In = atoms=='In'
def_P = atoms=='P'
def_Ga = atoms=="Ga"
def_InGa = np.logical_or(def_In,def_Ga)
def_cat=np.logical_or(def_InGa,atoms=="Al")
def_lig = np.logical_not(np.logical_or(def_cat,def_P))

all_def_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(coords,def_cat,def_P,def_lig)

#check for overlapping atoms
# for i,atom1 in enumerate(all_def_dists):
# 	double=False
# 	for j,dist in enumerate(atom1):
# 		if dist==0 and double==False:
# 			double=True
# 		elif dist==0 and double==True:
# 			print("we have an overlapping atom at overall index", i, "which is a", atoms[i])

#check that all hydrogens are in methyl groups
#for i,atom1 in enumerate(all_def_dists):
#	if atoms[i]=="H":
#		neighbor_count=0
#		for j,dist in enumerate(atom1):
#			if atoms[j]=="H" and dist<2:
#				neighbor_count=neighbor_count+1
#		if neighbor_count!=3:
#			print("Hyrdogen", to_atom_specific_index(atoms,i), "is not in a methyl group and instead has",neighbor_count, "neighbors")

#check for C not bound to another C
for i,atom1 in enumerate(all_def_dists):
	if atoms[i]=="C":
		c_bond=False
		for j,dist in enumerate(atom1):
			if dist<2 and atoms[j]=="C" and i!=j:
				c_bond=True
		if c_bond==False:
			print("officer the problem is carbon", to_atom_specific_index(atoms,i))



