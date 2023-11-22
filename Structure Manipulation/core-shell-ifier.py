import numpy as np
from matplotlib import pyplot as plt
from geom_helper import *
import sys
import math
import random

#the goal of this script is to be able to turn a core-only (unoptimized) QD into a core-shell QD of the same overall size
#the script substitutes Indium/Gallium for In/Ga/Al in the core and shell
#we want to be able to specify a number of monolayers for the shell

xyz=sys.argv[1] #the original unoptimized dot. InPF, InPCl, GaPF, GaPCl only
out=sys.argv[2] #name of .xyz to write
core_mat=sys.argv[3] #the desired core material, i.e. "InGaP" or "ingap" for InGaP. 
shell_mat=sys.argv[4] #the desired shell material, i.e. "InAlP" or "inalp" for InAlP
n_mono=int(sys.argv[5]) #the number of desired shell monolayers

cutoff = 3.05 #manual parameter for nearest neighbor bond distance

#note that core/shell materials have been manually implemented and some may not be supported


coords,atoms=read_input_xyz(xyz)

ind_in = atoms=='In'
ind_ga=atoms=="Ga"
ind_p = atoms=='P'
ind_F= atoms == 'F'
ind_Cl = atoms=="Cl"
ind_inga=np.logical_or(ind_in, ind_ga)
ind_attach=np.logical_or(ind_Cl,ind_F) #only F and Cl ligands supported for now
ind_lig=ind_attach
if np.count_nonzero(ind_p)==0:
	print("Current XYZ Anion Not Yet Supported")

dists, inp_dists, garbage, garbage2,pin_dists =get_dists(coords, ind_inga, ind_p)
in_underc_ind,p_underc_ind = get_underc_index(coords,ind_inga,ind_p,ind_lig,ind_attach,cutoff,4,verbose=False) #this gets us the indices of all P-3c, which is useful

connectivity=connectivity_finder(dists,cutoff)

surface=[] #determine which atoms are in the "outermost" shell. 1st loop just gets undercoordinated atoms and atoms bound to F
p_count=0
in_count=0
for i,atom in enumerate(atoms):
	if ind_p[i]:
		p_count=p_count+1
	if ind_inga[i]:
		in_count=in_count+1

	if atom=="P" and p_underc_ind[p_count-1]:
		surface.append(i)
	elif ind_inga[i] and in_underc_ind[in_count-1]:
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

	shell_cat=np.logical_and(shell_ind,ind_inga)
	shell_an=np.logical_and(shell_ind,ind_p)
	n_shell_cat=np.count_nonzero(shell_cat)
	n_shell_an=np.count_nonzero(shell_an)
	core_ind=np.logical_not(shell_ind)
	core_cat=np.logical_and(core_ind,ind_inga)
	core_an=np.logical_and(core_ind,ind_p)
	n_core_cat=np.count_nonzero(core_cat)
	n_core_an=np.count_nonzero(core_an)
	current_cat=atoms[np.nonzero(ind_inga)[0][0]] #what cation do we have now?

	print(n_shell_cat, n_shell_an, n_core_cat, n_core_an)

	if core_mat == "ingap" or core_mat=="InGaP" or core_mat=="gainp" or core_mat=="GaInP":		
		if current_cat=="Ga":
			replace_cat="In"
		elif current_cat=="In":
			replace_cat="Ga"
		else:
			print("Starting XYZ Core Cation Not Yet Supported")

		nsub=math.floor(n_core_cat/2)
		sub_ind=[]
		while len(sub_ind)<nsub: #generate random indices to substitute
			new_ind=random.randint(1,n_core_cat)
			if new_ind not in sub_ind:
				sub_ind.append(new_ind)

		print("Substituting", nsub, "out of", n_core_cat,"core", current_cat, "with", replace_cat)

		sub_count=0
		new_atoms=[]
		for i,atom in enumerate(atoms): #do the substitution
			if core_cat[i]:
				sub_count=sub_count+1
				if sub_count in sub_ind:
					atoms[i]=replace_cat


	if shell_mat=="inalp" or shell_mat=="InAlP" or shell_mat=="alinp" or shell_mat=="AlInP":
		replace_cat1="In"
		replace_cat2="Al"

		nsub=math.floor(n_shell_cat/2)
		sub_ind=[]
		while len(sub_ind)<nsub: #generate random indices to substitute
			new_ind=random.randint(1,n_shell_cat)
			if new_ind not in sub_ind:
				sub_ind.append(new_ind)

		print("Substituting", nsub, "out of", n_shell_cat,"shell", current_cat, "with", replace_cat1, "and the rest with", replace_cat2)

		sub_count=0
		new_atoms=[]
		for i,atom in enumerate(atoms): #do the substitution
			if shell_cat[i]:
				sub_count=sub_count+1
				if sub_count in sub_ind:
					atoms[i]=replace_cat1
				else:
					atoms[i]=replace_cat2


	if shell_mat=="znse" or shell_mat=="ZnSe" or shell_mat=="ZnS" or shell_mat=="zns":
		replace_cat="Zn"
		if len(shell_mat)==3:
			replace_an="S"
		else:
			replace_an="Se"

		for i,atom in enumerate(atoms):
			if shell_cat[i]:
				atoms[i]=replace_cat
			if shell_an[i]:
				atoms[i]=replace_an

		#still need to do some ligand finagling to reach an acceptable charge
		charge=0
		for i,atom in enumerate(atoms):
			if atom=="In" or atom=="Ga" or atom =="Al":
				charge = charge+3
			elif atom =="P":
				charge = charge-3
			elif atom=="F" or atom =="Cl":
				charge = charge-1
			elif atom=="Zn":
				charge = charge+2
			elif atom=="Se" or atom=="S":
				charge = charge-2
			else:
				print("You have an element in here that isn't supported!")
		print("Charge is ", charge, " after ZnSe substitution")

write_xyz(out,atoms,coords)



#for testing
#new_atoms=[]
#new_coords=[]
#for i,atom in enumerate(atoms):
#	if i not in shell:
#		new_atoms.append(atom)
#		new_coords.append(coords[i])

#write_xyz(out, new_atoms, new_coords)
