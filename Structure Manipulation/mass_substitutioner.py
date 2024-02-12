import numpy as np
import sys
import math
import random
from geom_helper import *
import os

# similar to random_substitutioner, writes many substituted versions of a base .xyz to subdirectories
# takes an .xyz and makes X copies of it, each with a unique Y of element A substituted to element B
# doesn't try to adjust bond lengths or anything

#inputs
xyz=sys.argv[1] #.xyz with original structure. assumed to be of the form X_DESCRIPTOR.xyz
orig=sys.argv[2] # element to replace
new=sys.argv[3] #element to add
nstruct=int(sys.argv[4]) #number of structures to generate
nsub=float(sys.argv[5]) #number of substitutions. If less than 1, treated as a fraction

coords,atoms=read_input_xyz(xyz)

flag=0
descriptor=""
for i,char in enumerate(xyz):
	if char=="_" and flag==0:
		flag=1
	elif char==".":
		flag=0
	elif flag==1:
		descriptor=descriptor+char

print(descriptor)

structure_count=0
substitution_list=[] #ensures each random copy is different
while structure_count<nstruct:

	#first, choose a set of nsub random orig atoms to substitute
	ind_orig = atoms==orig
	n_orig=np.count_nonzero(ind_orig)

	if nsub<1: #in case nsub is a fraction
		nsub=int(round(nsub*n_orig))
	else:
		nsub=int(nsub)

	sub_candidate=random.sample(range(n_orig),nsub) #random set of indices, without repeats
	sub_candidate.sort()

	if sub_candidate not in substitution_list: #ensures uniqueness
		substitution_list.append(sub_candidate)
		sub_candidate=np.array(sub_candidate)

		#the for loop way

		# new_atoms=[]
		# orig_count=0
		# for i,atom in enumerate(atoms):
		# 	if atom==orig:
		# 		if orig_count in sub_candidate:
		# 			new_atoms.append(new)
		# 		else:
		# 			new_atoms.append(atom)
		# 		orig_count=orig_count+1
		# 	else:
		# 		new_atoms.append(atom)

		#faster way
		overall_sub_indices = np.where(ind_orig)[0][sub_candidate]
		new_atoms=np.copy(atoms)
		new_atoms[overall_sub_indices]=new

		#write file
		structure_count=structure_count+1
		dir_name=str(structure_count)+"_sub"
		if not os.path.exists(dir_name):
			os.mkdir(dir_name)

		write_xyz(dir_name+"/unopt_"+dir_name+"_"+descriptor+".xyz",new_atoms,coords)

		





		


	



