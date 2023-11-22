import sys
from geom_helper import *

#a quick script to print the atom-specific indexes of all species bound to a given atomic type
# to be used in doped systems

xyz = sys.argv[1] #the .xyz file
dopant=sys.argv[2] #the dopant species you want indices for

threshold=3 #somewhat arbitrary, to determine if two atoms are bonded. May need to be changed for certain systems

coords,atoms=read_input_xyz(xyz)

#get all bond distances
#these don't really matter
ind_In = atoms=='In'
ind_P = atoms=='P'
ind_Ga = atoms=="Ga"
ind_InGa = np.logical_or(ind_In,ind_Ga)
ind_lig = np.logical_not(np.logical_or(ind_InGa,ind_P))
lig="Cl" 

all_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(coords,ind_InGa,ind_P,ind_lig) #don't care about other dists, just an artefact of the helper I'm using

connectivity=connectivity_finder(all_dists,threshold)

targets=[]
for i,atom in enumerate(atoms): 
	if atom==dopant:
		targets.extend(connectivity[i])

targets=list(set(targets)) #removes duplicates

species=[]
for i,target in enumerate(targets):
	if atoms[target] not in species:
		species.append(atoms[target])

for i,element in enumerate(species):
	print(element,end=" ")
	for i,target in enumerate(targets):
		if atoms[target]==element:
			print(to_atom_specific_index(atoms,target),end=" ")
	print()


