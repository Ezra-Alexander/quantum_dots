import numpy as np
import sys
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *
import random

#the brute's way of passivating a QD from a negative charge to neutral
#goes through and removes a certain number of -1 halogen ligands from 4c surface cations
#i have updated this to avoid the creation of strong dipole moments

xyz=sys.argv[1] #the original dot. In, P, Ga, Al, F, Cl supported
out=sys.argv[2] #name of .xyz to write
charge=int(sys.argv[3]) #the total negative charge of the dot, i.e. the total number of ligands to remove

cutoff = 2.8 #manual parameter for cation-ligand distance. May need to be tweaked
nncutoff=3 #We want to avoid making any In-2c

coords,atoms = read_input_xyz(xyz)

n_pop=0
while n_pop < charge:	

	#ind_p = np.logical_or(atoms=='P',np.logical_or(atoms=="Se",atoms=="S"))
	ind_p = np.logical_or(np.logical_or(atoms=='P',np.logical_or(atoms=="Se",atoms=="S")),np.logical_or(atoms=="C",atoms=="Si")) #temporary for P to C/Si defects
	ind_cat = np.logical_or(atoms=="In",np.logical_or(atoms=="Ga", np.logical_or(atoms=="Al",atoms=="Zn")))
	ind_lig= np.logical_or(atoms=="F", atoms=="Cl")
	ind_attach = ind_lig #relic

	#print(ind_cat)

	#what's undercoordinated?
	cat_underc_ind,p_underc_ind = get_underc_index(coords,ind_cat,ind_p,ind_lig,ind_attach,cutoff,nncutoff,verbose=False)
	n_underc_cat=np.count_nonzero(cat_underc_ind) #the number of in-2c
	#print(n_underc_cat)

	n_lig=np.count_nonzero(ind_lig)

	target=random.randint(1,n_lig)

	count=0
	for i,atom in enumerate(atoms):
		if ind_lig[i]:
			count=count+1
			if count==target:
				target_total_index=i


	#add a second target that is maximal distance from the first target
	if (charge-n_pop)>1:
		dists=dist_all_points(coords)
		connectivity=connectivity_finder(dists,cutoff)
		max_dist=0
		count=0
		for i,atom in enumerate(atoms):
			if ind_lig[i]:
				count=count+1
				if dists[i][target_total_index]>max_dist:
					anchor_i=connectivity[i][0]
					if len(connectivity[anchor_i])>3:
						second_target=count
						max_dist=dists[i][target_total_index]
	print(target, second_target)

	new_coords=[]
	new_atoms=[]
	lig_count=0
	for j,atom in enumerate(atoms):
		if ind_lig[j]:
			lig_count=lig_count+1
			if lig_count==target or lig_count==second_target:
				pass
			else:
				new_coords.append(coords[j])
				new_atoms.append(atoms[j])
		else:
			new_coords.append(coords[j])
			new_atoms.append(atoms[j])

	new_coords=np.array(new_coords)
	new_atoms=np.array(new_atoms)

	#ind_p = np.logical_or(new_atoms=='P',np.logical_or(new_atoms=="Se",new_atoms=="S"))
	ind_p = np.logical_or(np.logical_or(new_atoms=='P',np.logical_or(new_atoms=="Se",new_atoms=="S")),np.logical_or(new_atoms=="C",new_atoms=="Si")) #temporary for P to C/Si defects
	ind_cat = np.logical_or(new_atoms=="In",np.logical_or(new_atoms=="Ga", np.logical_or(new_atoms=="Al",new_atoms=="Zn")))
	ind_lig= np.logical_or(new_atoms=="F", new_atoms=="Cl")
	ind_attach = ind_lig #relic

	new_cat_underc_ind,new_p_underc_ind = get_underc_index(new_coords,ind_cat,ind_p,ind_lig,ind_attach,cutoff,nncutoff,verbose=False)
	new_n_underc_cat=np.count_nonzero(new_cat_underc_ind)
	#print(new_n_underc_cat)

	if new_n_underc_cat==n_underc_cat:
		coords=new_coords
		atoms=new_atoms
		if (charge-n_pop)>1:
			n_pop=n_pop+2
		else:
			n_pop=n_pop+1
		print("Removed", n_pop, "ligands")




write_xyz(out,atoms,coords)




