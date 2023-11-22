import numpy as np
import sys
import math
from qchem_helper import read_xyz, write_xyz
from geom_helper import *

#This code is meant to compute an "average inverse P to ligand distance" for P in InP QDs in order to rationalize XPS Chemical Shifts
#Didn't end up being helpful

#The orginal .xyz file that I am editing
xyz = sys.argv[1]

#The phosphorous index to compute distances for
phos = sys.argv[2]

coords, atoms = read_xyz(xyz)

ind_In = atoms=='In'
ind_P = atoms=='P'
ind_Ga = atoms=="Ga"
ind_InGa = np.logical_or(ind_In,ind_Ga)
ind_InP= np.logical_or(ind_In,ind_P)
ind_lig = (atoms == "Br")                     #change this line when the ligand changes

all_dists,inp_dists,inf_dists,inpf_dists,pin_dists = get_dists(coords,ind_In,ind_P,ind_lig)

index = 0
pcount = 0
for i,atom in enumerate(atoms):
	if atom == "P":
		pcount=pcount+1
		if pcount==int(phos):
			break
	index = index+1

print(index)

inv_dist = 0
for i,dist in enumerate(all_dists[index]):
	if ind_lig[i] and dist <6:
		inv_dist = inv_dist + (1/dist)

print(inv_dist)