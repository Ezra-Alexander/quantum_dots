import numpy as np
from matplotlib import pyplot as plt
from geom_helper import *
import sys
import math
import random

#takes a qd .xyz and prints out info about its core-shell makeup, size, etc.

xyz=sys.argv[1] #xyz of dot
n_mono=int(sys.argv[2]) #of monolayers of shell this dot is expected to have, generous definition

cutoff = 2.9 #manual parameter for nearest neighbor bond distance

coords,atoms=read_input_xyz(xyz)

ind_in = atoms=='In'
ind_ga=atoms=="Ga"
ind_al=atoms=="Al"
#ind_al=atoms=="Zn"
ind_p = atoms=='P'
ind_p=np.logical_or(atoms=="P",atoms=="S")
ind_F= atoms == 'F'
ind_Cl = atoms=="Cl"
ind_inga=np.logical_or(ind_in, ind_ga)
ind_cat=np.logical_or(ind_inga, ind_al)
ind_attach=np.logical_or(ind_Cl,ind_F) #only F and Cl ligands supported for now
ind_lig=ind_attach
if np.count_nonzero(ind_p)==0:
	print("Current XYZ Anion Not Yet Supported")

dists, inp_dists, garbage, garbage2,pin_dists =get_dists(coords, ind_cat, ind_p)
in_underc_ind,p_underc_ind = get_underc_index(coords,ind_cat,ind_p,ind_lig,ind_attach,cutoff,4,verbose=False) #this gets us the indices of all P-3c, which is useful

connectivity=connectivity_finder(dists,cutoff)

surface=[] #determine which atoms are in the "outermost" shell. 1st loop just gets undercoordinated atoms and atoms bound to F
p_count=0
in_count=0
for i,atom in enumerate(atoms):
	if ind_p[i]:
		p_count=p_count+1
	if ind_cat[i]:
		in_count=in_count+1

	if atom=="P" and p_underc_ind[p_count-1]:
		surface.append(i)
	elif ind_cat[i] and in_underc_ind[in_count-1]:
		surface.append(i)
	else:
		for j,bond in enumerate(connectivity[i]):
			if ind_lig[bond]:
				surface.append(i)


#now the 2nd loop gets everything bound to the "surface"
shell=[]
for i,atom in enumerate(atoms):
	if i in surface:
		shell.append(i)
	else:
		for j,bond in enumerate(connectivity[i]):
			if bond in surface and (i not in shell):
				shell.append(i)

while n_mono > 1: #do subsequent monolayers of shell
	shell2=[]
	for i,atom in enumerate(atoms):
		if i in shell:
			shell2.append(i)
		else:
			for j,bond in enumerate(connectivity[i]):
				if bond in shell and (i not in shell2):
					shell2.append(i)

	#shell3=[]
	#for i,atom in enumerate(atoms):
	#	if i in shell2:
	#		shell3.append(i)
	#	else:
	#		for j,bond in enumerate(connectivity[i]):
	#			if bond in shell2 and (i not in shell3):
	#				shell3.append(i)

	#shell=shell3
	shell=shell2

	n_mono=n_mono-1

if n_mono <= 1: #if n_mono = 1, we are done selecting our "shell"
	shell_ind=[]
	for i,atom in enumerate(atoms): #get shell indexes
		if i in shell:
			shell_ind.append(True)
		else:
			shell_ind.append(False)
	np.array(shell_ind)

	shell_cat=np.logical_and(shell_ind,ind_cat)
	shell_an=np.logical_and(shell_ind,ind_p)
	shell_in=np.logical_and(shell_ind,ind_in)
	shell_ga=np.logical_and(shell_ind,ind_ga)
	shell_al=np.logical_and(shell_ind,ind_al)
	n_shell_cat=np.count_nonzero(shell_cat)
	n_shell_in=np.count_nonzero(shell_in)
	n_shell_ga=np.count_nonzero(shell_ga)
	n_shell_al=np.count_nonzero(shell_al)
	n_shell_an=np.count_nonzero(shell_an)

	core_ind=np.logical_not(shell_ind)
	core_cat=np.logical_and(core_ind,ind_cat)
	core_in=np.logical_and(core_ind,ind_in)
	core_ga=np.logical_and(core_ind,ind_ga)
	core_al=np.logical_and(core_ind,ind_al)
	core_an=np.logical_and(core_ind,ind_p)
	n_core_cat=np.count_nonzero(core_cat)
	n_core_in=np.count_nonzero(core_in)
	n_core_ga=np.count_nonzero(core_ga)
	n_core_al=np.count_nonzero(core_al)
	n_core_an=np.count_nonzero(core_an)

	print("The core consists of", n_core_in, 'In', n_core_ga, "Ga", n_core_al, "Al, and", n_core_an, "P")
	print("The shell consists of", n_shell_in, 'In', n_shell_ga, "Ga", n_shell_al, "Al", n_shell_an, "P, and", np.count_nonzero(ind_lig), "F")

	print("It has charge", (np.count_nonzero(ind_cat)*3)-(np.count_nonzero(ind_p)*3)-np.count_nonzero(ind_lig))
	
	np.save('in_core_ind',core_in)
	np.save('p_core_ind',core_an)

