import numpy as np
import sys
import math
from geom_helper import *
import random

#UNFINISHED

#for a given QD .xyz file, this replaces all the halogen ligands with oxygens, and then grows dangling carboxyl groups off of each of them
#the hope is that the geometry optimization will sort everything out
#ideally this should be run on an unoptimized starting structure

xyz=sys.argv[1] #the initial .xyz file
out_file=sys.argv[2] #what you want to name the output .xyz file

threshold=2.6 #manual parameter determining upper bound for ligand-cation bond length

coords,atoms = read_input_xyz(xyz)

ind_In = atoms=='In'
ind_P = atoms=='P'
ind_Ga = atoms=="Ga"
ind_InGa = np.logical_or(ind_In,ind_Ga)
ind_cat=np.logical_or(ind_InGa, np.logical_or(atoms=="Al", atoms=="Zn"))
ind_an=np.logical_or(ind_P, np.logical_or(atoms=="S",atoms=="Se"))
ind_lig = np.logical_not(np.logical_or(ind_cat,ind_an))

all_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(coords,ind_InGa,ind_P,ind_lig) #don't care about other dists, just an artefact of the helper I'm using

#manual carboxyl geometry, from magic cluster
#reference_atoms= np.array(["O", "O", "C", "C", "H", "H", "H"])
#reference_coords= np.array([[-6.8897180187,        1.0021509810,       3.6491191817], [-8.6171143274,        0.7187599606,       2.2710465014],[-8.1372668061,        0.8355114352,       3.4514910519],[-9.0984962138,        0.7858724048,       4.6143303048],[-8.7163149861,        0.0954854001,       5.3917381764],[-10.1120036672,        0.4989351368,       4.2848847320],[-9.1279297603,       1.7911254084,       5.0834503812]])

final_atoms=[]
final_coords=[]
for i,atom in enumerate(atoms):
	if ind_lig[i]:
		final_atoms.append("O") # do swap
		final_coords.append(coords[i])

		#get bonded cations
		coordination = np.count_nonzero(all_dists[i]<threshold) #note this will count itself
		bonds=np.where(all_dists[i]<threshold)[0]
		bonds=np.delete(bonds,np.where(bonds==i))

		if coordination==2: #terminal ligands
			o1=coords[i]
			cat=coords[bonds[0]]
			ref_bond=o1-cat
			





	else:
		final_atoms.append(atom)
		final_coords.append(coords[i])


write_xyz(out_file,final_atoms,final_coords)