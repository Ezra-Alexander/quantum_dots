import sys
from geom_helper import *

#instead of painstakingly entering every index we want, what if we could just write a set of indexes and get some number of nearest neighbors

xyz = sys.argv[1] 

num_nn=int(sys.argv[2])
nn_name=num_nn

num_list = [] #still overall indices
for num in sys.argv[3:]:
	num_list.append(int(num)-1) #0-centering

threshold=2.9 #somewhat arbitrary, to determine if two atoms are bonded. May need to be changed for certain systems

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


to_carve=num_list.copy()
while num_nn>0:

	temp=[]
	for i,atom in enumerate(to_carve):
		for j,bond in enumerate(connectivity[atom]):
			if (bond not in to_carve) and (bond not in temp):
				temp.append(bond)

	to_carve.extend(temp)
	num_nn=num_nn-1

new_coords=[]
new_atoms=[]
for i,atom in enumerate(atoms):
	if i in to_carve:
		new_atoms.append(atom)
		new_coords.append(coords[i])

write_xyz(str(nn_name)+"_pics.xyz", new_atoms, new_coords)

