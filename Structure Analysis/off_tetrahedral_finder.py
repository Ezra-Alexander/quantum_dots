import numpy as np
import sys
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *

#This script is meant to read a .xyz file and print out each non-ligand 4-coordinate atom that has a strongly off-tetrahedral geometry (1 bond angle >= 20 degrees from 109.5)
#not properly configured for Gallium yet
#also doesn't work for ligands with more than one atom (yet)

#This part all pretty much comes from geom_analysis_ezra.py and determines which atoms are < 4c (we'll include 5c atoms for now)
print()
cutoff = 2.85  # nearest neighbor cutoff distance (lowest) #USER SPECIFIED
print('Bond distance cutoff: ',cutoff)
nncutoff = 4  # number of nearest neighbors to be considered "passivated" (incl. ligands) #USER SPECIFIED
angle_cutoff=8 #minimum degrees from 109.5 for an angle to be considered "off-tetrahedral" #USER SPECIFIED
print('Angle cutoff: ',angle_cutoff)
pp_cutoff=4.1 #cutoff for two P to be "adjacent" in Angstroms #USER SPECIFIED
print('Adjacent P cutoff: ',pp_cutoff)
lig_atom = "O" #USER SPECIFIED
print('Ligand: ',lig_atom)
verbose = 'True'==sys.argv[2] #supress of allow printing of each off-tet In and P

QD_file=sys.argv[1]

QD_xyz,atom_names = read_input_xyz(QD_file)

ind_In = atom_names=='In'
ind_P = atom_names=='P'
ind_Ga = atom_names=="Ga"
ind_InGa = np.logical_or(ind_In,ind_Ga)
ind_InP= np.logical_or(ind_In,ind_P)
ind_lig = (atom_names == lig_atom)
ind_attach = (atom_names == lig_atom)

print(np.count_nonzero(ind_In),"In", np.count_nonzero(ind_P),"P", np.count_nonzero(ind_Ga),"Ga")
print()

in_underc_ind_s,p_underc_ind_1 = get_underc_index(QD_xyz,ind_In,ind_P,ind_lig,ind_attach,cutoff,nncutoff,verbose=False)
if np.count_nonzero(ind_Ga) > 0:
    ga_underc_ind_s, pee = get_underc_index(QD_xyz,ind_Ga,ind_P,ind_lig,ind_attach,cutoff,nncutoff,verbose=False)
    poo, p_underc_ind_2 = get_underc_index(QD_xyz,ind_InGa,ind_P,ind_lig,ind_attach,cutoff,nncutoff,verbose=False)


#print('Starting geometry')
print('Undercoordinated In:',np.count_nonzero(in_underc_ind_s))
for i,boul in enumerate(in_underc_ind_s):
	if boul:
		print("Indium ", i+1)
print()
print('Undercoordinated P (no Ga):',np.count_nonzero(p_underc_ind_1))
for i,boul in enumerate(p_underc_ind_1):
	if boul:
		print("Phosphorous ", i+1)
print()

if np.count_nonzero(ind_Ga) > 0:
    print('Undercoordinated Ga:',np.count_nonzero(ga_underc_ind_s))
    print('Undercoordinated P w/ Ga:',np.count_nonzero(p_underc_ind_2))
    print()


#Do the thing!
dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(QD_xyz,ind_In, ind_P,ind_lig)
in_off_tet_ind, p_off_tet_ind, in_angles, p_angles = get_off_tet_index(QD_xyz, ind_In, ind_P,ind_lig,in_underc_ind_s,p_underc_ind_1,dists,cutoff,angle_cutoff)
#for Ga #in_off_tet_ind, p_off_tet_ind, in_angles, p_angles = get_off_tet_index(QD_xyz, ind_In, ind_P,ind_lig,in_underc_ind_s,p_underc_ind_2,dists,cutoff,angle_cutoff)

print("This structure has a total of", sum(in_off_tet_ind), "off-tetrahedral In and",sum(p_off_tet_ind),"off-tetrahedral P.")
print()

if verbose:
	print("The following 4-Coorinate Indium have at least 1 off-tetrahedral bond angle:")
	for i,indium in enumerate(in_off_tet_ind):
		if indium:
			print()
			print("Indium",str(i+1) + ":")
			print("Which has the following bond angles:")
			for j,angle in enumerate(in_angles[i]):
				print(str(angle))

	print()
	print("The following 4-Coorinate Phosphorous have at least 1 off-tetrahedral bond angle:")
	for i,phos in enumerate(p_off_tet_ind):
		if phos:
			print()
			print("Phosphorous",str(i+1) + ":")
			print("Which has the following bond angles:")
			for j,angle in enumerate(p_angles[i]):
				print(str(angle))

p_pairs = get_off_pairs(dists,p_off_tet_ind,pp_cutoff,ind_P)

print(p_pairs)



