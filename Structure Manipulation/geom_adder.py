import numpy as np
import sys
import math
from geom_helper import *
from qchem_helper import read_xyz, write_xyz

#The orginal .xyz file that I am editing
xyz = sys.argv[1]

#The name of the .xyz file this will make
name = sys.argv[2]

#The overall index of the 3-coordinate atom to add something to
index = int(sys.argv[3])-1

#what atom I want to add
add = sys.argv[4]

coords, atoms = read_xyz(xyz)

#find vector along 4th dimension
ind_In = atoms=='In'
ind_P = atoms=='P'
ind_Ga = atoms=="Ga"
ind_InGa = np.logical_or(ind_In,ind_Ga)
ind_cat=np.logical_or(ind_InGa,atoms=="Al")
ind_InP= np.logical_or(ind_In,ind_P)
ind_lig = (atoms == "Cl")

all_dists,inp_dists,inf_dists,inpf_dists,pin_dists = get_dists(coords,ind_cat,ind_P,ind_lig)
#all_dists,inp_dists,inf_dists,inpf_dists,pin_dists = get_dists(coords,ind_Ga,ind_P,ind_lig)

target_dists = all_dists[index]
#sorted_dists = np.sort(target_dists)
indexes = [1000,1000,1000]
dists = [1000,1000,1000]
for i,dist in enumerate(target_dists):
	if dist < max(dists) and i!=index:
		indexes.pop(dists.index(max(dists)))
		dists.pop(dists.index(max(dists)))
		indexes.append(i)
		dists.append(dist)

ind_1 = indexes[0]
ind_2 = indexes[1]
ind_3 = indexes[2]

#bond_1 = sorted_dists[1]
#bond_2 = sorted_dists[2]
#bond_3 = sorted_dists[3]

#this assumes that the three closest dists are not exactly equal
#ind_1 = np.where(target_dists==bond_1)[0][0]
#ind_2 = np.where(target_dists==bond_2)[0][0]
#ind_3 = np.where(target_dists==bond_3)[0][0]

dummy = (coords[ind_1] + coords[ind_2] + coords[ind_3])/3

#dummy = (coords[ind_2] + coords[ind_3])/3 #temp

vector = coords[index]-dummy
mag = np.linalg.norm(vector)
unit = vector/mag

#Need to add bond lengths of each type as I go
if add == "In" and atoms[index] == "P":
	new_coords = coords[index] + unit*2.58
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add == "F" and atoms[index] == "In":
	new_coords = coords[index] + unit*2.2
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add == "O" and atoms[index] == "In":
	new_coords = coords[index] + unit*2.05
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add == "Ga" and atoms[index] == "P":
	new_coords = coords[index] + unit*2.36
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add=="F" and atoms[index]=="Ga":
	new_coords = coords[index] + unit*2.12
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add=="Al" and atoms[index]=="P":
	new_coords = coords[index] + unit*2.32
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add=="Cl" and atoms[index]=="In":
	new_coords = coords[index] + unit*2.36
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add=="Cl" and atoms[index]=="Zn":
	new_coords = coords[index] + unit*2.186
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)
if add=="Cl" and atoms[index]=="Al":
	new_coords = coords[index] + unit*2.116
	print("Done!")
	new_atoms = np.append(atoms,add)
	final_coords = np.append(coords,[new_coords], axis=0)


write_xyz(name, new_atoms, final_coords)

