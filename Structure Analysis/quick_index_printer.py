import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *

#quickly print some indices I can use to run average_param_getter
#let's just hard-code 4c for now

xyz=sys.argv[1] #the xyz file in question
element=sys.argv[2] # just In or P so far
underc=int(sys.argv[3]) #do you want all the underc or all the 4c? enter the coordination # you want. only does 3 or 4

cutoff = 2.9  # nearest neighbor cutoff distance (lowest)
lig_atom = "F"
nncutoff = 4 

#info the code has
masses = {"Ga":69.723, "P":30.974, "In":114.82,"Cl":35.453,"F":18.998403162}

coords,atoms = read_input_xyz(xyz)

ind_Cd = atoms=='In'
ind_Se = atoms=='P'
ind_Ga = atoms=="Ga"
ind_Al = atoms=="Al"    
ind_InGa = np.logical_or(ind_Cd,ind_Ga)
ind_cat=np.logical_or(ind_InGa,ind_Al)
ind_CdSe= np.logical_or(ind_Cd,ind_Se)
ind_F= (atoms == lig_atom)
ind_O = atoms=="O"
ind_attach=np.logical_or(ind_O,ind_F)
ind_lig=ind_attach

in_underc, p_underc=get_underc_index(coords,ind_cat,ind_Se,ind_lig,ind_attach,cutoff,nncutoff,verbose=False)

if underc==4:
	if element=="In":
		in_count=0
		for i,atom in enumerate(atoms):
			if ind_Cd[i]:
				in_count=in_count+1
				if not in_underc[in_count-1]:
					print(in_count,end=" ")
	elif element=="P":
		p_count=0
		for i,atom in enumerate(atoms):
			if ind_Se[i]:
				p_count=p_count+1
				if not p_underc[p_count-1]:
					print(p_count, end=" ")
	else:
		print("Requested element not yet supported")
if underc<4:
	if element=="In":
		in_count=0
		for i,atom in enumerate(atoms):
			if ind_Cd[i]:
				in_count=in_count+1
				if in_underc[in_count-1]:
					print(in_count,end=" ")
	elif element=="P":
		p_count=0
		for i,atom in enumerate(atoms):
			if ind_Se[i]:
				p_count=p_count+1
				if p_underc[p_count-1]:
					print(p_count, end=" ")
	else:
		print("Requested element not yet supported")
