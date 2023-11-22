import numpy as np
import sys
import math
from geom_helper import *


#like the general carver, but now designed to work for 3c In and Ga
#I want this to specifically work for one .xyz, and two In-3c in that .xyz
#and this is a single atom interpolation, InF3 or GaF3 only
#which is to say that we will set all bonded P to be F and adjust the bond length accordingly


xyz=sys.argv[1] #.xyz
target_species=sys.argv[2] 
start_index=int(sys.argv[3]) #atom specific, any number
end_index=int(sys.argv[4])

threshold=2.9 #somewhat arbitrary, to determine if two atoms are bonded. May need to be changed for certain systems

if start_index==end_index:
	print("I guess you can do an interpolation from something to itself, but I don't think it will be very interesting ...")


#read xyzs
coords, atoms = read_input_xyz(xyz)


#extract the targeted atoms

#convert to overall indices
start_target=0
end_target=0
count = 0
for i,atom in enumerate(atoms):
	if atom==target_species:
		count=count+1
		if count==start_index:
			start_target=i
		if count==end_index:
			end_target=i


#get all bond distances
ind_In = atoms=='In'
ind_P = atoms=='P'
ind_Ga = atoms=="Ga"
ind_InGa = np.logical_or(ind_In,ind_Ga)
ind_lig = np.logical_not(np.logical_or(ind_InGa,ind_P))
lig="F" 

all_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(coords,ind_InGa,ind_P,ind_lig) #don't care about other dists, just an artefact of the helper I'm using


#now pull out everything bonded to your targets
bonded_start=[]
bonded_start_atom_order=[]
bonded_end=[]
bonded_end_atom_order=[]
for i,atom in enumerate(atoms):
	if atom!=target_species:
			if all_dists[i][start_target]<threshold:
				if i not in bonded_start:
					bonded_start.append(i)
					bonded_start_atom_order.append(atom)

			if all_dists[i][end_target]<threshold:
				if i not in bonded_end:
					bonded_end.append(i)
					bonded_end_atom_order.append(atom)


#check to make sure both are three coordinate
if len(bonded_start)!= len(bonded_end):
	print("Problem! The coordination of your two selected atoms is different")
	print("You start with",len(bonded_start),"and end with", len(bonded_end))


#now I substitute all P to F

#some futzing around with bond lengths
f_bond_length=1.986 #terminal In-F bond length
ideal_bond=2.579 #from bulk crystal
if target_species=="Ga":
		f_bond_length=1.774 #terminal Ga-F bond length
		ideal_bond=2.360 #from bulk crystal
if lig!="F":
	print("Error - You need to hard code this cation-ligand bond length")

start_coords=np.copy(coords)
end_coords=np.copy(coords)


#for the start geometry
for i,atom in enumerate(bonded_start_atom_order):
	if atom=="P":
		atoms[bonded_start[i]]=lig
		bond=coords[bonded_start[i]]-coords[start_target]
		bond_length=np.linalg.norm(bond)
		normal=bond/bond_length
		temp_coords=coords[start_target]+(normal*f_bond_length*(bond_length/ideal_bond))
		start_coords[bonded_start[i]]=temp_coords

#for the final geometry
for i,atom in enumerate(bonded_end_atom_order):
	if atom=="P":
		atoms[bonded_end[i]]=lig
		bond=coords[bonded_end[i]]-coords[end_target]
		bond_length=np.linalg.norm(bond)
		normal=bond/bond_length
		temp_coords=coords[end_target]+(normal*f_bond_length*(bond_length/ideal_bond))
		end_coords[bonded_end[i]]=temp_coords		

#one remaining problem
#we need to rotate everything in space so that the two structures are as close to each other as possible

#here we switch the indices in the final structure so that the F are as close as possible to their corresponding intial structure counterparts
for j,start_atom in enumerate(bonded_start):
	current_dist=np.linalg.norm(end_coords[bonded_end[j]]-start_coords[start_atom]-end_coords[end_target]+start_coords[start_target])
	min_dist=1000
	min_i=10	
	for i,end_atom in enumerate(bonded_end):
		distance=np.linalg.norm(end_coords[end_atom]-start_coords[start_atom]-end_coords[end_target]+start_coords[start_target])
		if i!=j and distance<min_dist:
			min_dist=distance
			min_i=i
	if j!=min_i and min_dist<current_dist: 
		end_coords[[bonded_end[j],bonded_end[min_i]]]=end_coords[[bonded_end[min_i], bonded_end[j]]]



			
n_atoms=1+len(bonded_start)
#now, write endpoints.xyz
with open("endpoints.xyz","w") as out:
	out.write(str(n_atoms)+"\n")
	out.write("Start \n")

	#starting center cation
	out.write(target_species)
	out.write("  ")
	out.write(str(start_coords[start_target][0]))
	out.write("  ")
	out.write(str(start_coords[start_target][1]))
	out.write("  ")
	out.write(str(start_coords[start_target][2]))
	out.write(" \n")

	#BONDED TO starting CENTER CATION
	for i,index in enumerate(bonded_start):
			out.write(str(atoms[index]))
			out.write("  ")
			out.write(str(start_coords[index][0]))
			out.write("  ")
			out.write(str(start_coords[index][1]))
			out.write("  ")
			out.write(str(start_coords[index][2]))
			out.write(" \n")


	#final cation
	out.write(str(n_atoms)+"\n")
	out.write("End \n")

	out.write(target_species)
	out.write("  ")
	out.write(str(end_coords[end_target][0]))
	out.write("  ")
	out.write(str(end_coords[end_target][1]))
	out.write("  ")
	out.write(str(end_coords[end_target][2]))
	out.write(" \n")

	for i,index in enumerate(bonded_end):
		out.write(str(atoms[index]))
		out.write("  ")
		out.write(str(end_coords[index][0]))
		out.write("  ")
		out.write(str(end_coords[index][1]))
		out.write("  ")
		out.write(str(end_coords[index][2]))
		out.write(" \n")
		