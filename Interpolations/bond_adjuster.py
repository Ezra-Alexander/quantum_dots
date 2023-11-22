import numpy as np
import sys
import math
from geom_helper import *

#this code is meant to take all bonds of a given type (i.e P-In bonds) in a given structure and adjust them to a new bond length
#as written this only works for two element systems
#doesn't work with 4-membered rings
#also note that it doesn't actually change any atoms to a different element

#need to adjust this in a few ways 
#so that it works with double .xyzs
#so that it actually changes the majority element to Lithium

xyz = sys.argv[1] #original .xyz file

name = sys.argv[2] #of the new xyz file

element1 = sys.argv[3] #the test is coded such that the minority element is element 2 for some reason
element2 = sys.argv[4]

new_length = sys.argv[5]

cutoff = 2.9 #set manually, may need to be changed. #make sure it is greater than the minimum distance between identical elements

coords,atoms = read_input_xyz(xyz)

ind_1 = atoms==element1
ind_2 = atoms==element2

dists, a, b, c, d = get_dists(coords, ind_1, ind_2)

elem_1_angles, elem_2_angles = get_angles(coords, ind_1, ind_2, dists, cutoff)

main_elem = "H" #this is a little confusing - the "main" element is just the one that's more prevalent, and the "min" element is often more important
min_elem = "C" #these settings also don't mean anything - I'm just initializing
if np.count_nonzero(ind_1) > np.count_nonzero(ind_2):
	main_elem=element1
	min_elem=element2
else:
	main_elem=element2
	min_elem=element1

new_coords = [0]
new_atoms=[]
transfer_atoms=[]
done_atoms=[]
count = 0
coords2 = coords.copy()

for i, atom in enumerate(atoms):
	if atom == min_elem and count==0:
		new_coords[0]=list(coords[i])
		new_atoms.append(min_elem)
		count = count+1
		done_atoms.append(i)

		for j,atom2 in enumerate(atoms):
			if atom2==main_elem and dists[i][j] < cutoff:
				vector = coords[j]-coords[i]
				unit = vector/np.linalg.norm(vector)
				coords_new = coords[i]+unit*float(new_length)
				new_coords.append(list(coords_new))
				coords2[j]=coords_new
				new_atoms.append(main_elem)
				transfer_atoms.append(j)
				done_atoms.append(j)

	elif atom == min_elem and count>0:
		
		count=count+1
		for j,atom2 in enumerate(atoms):
			if atom2==main_elem and dists[i][j] < cutoff and j in transfer_atoms:				
				vector = coords[i]-coords[j]
				unit = vector/np.linalg.norm(vector)
				coords_new = coords2[j]+unit*float(new_length)
				new_coords.append(list(coords_new))
				coords2[i]=coords_new
				new_atoms.append(min_elem)
				done_atoms.append(i)

		transfer_atoms=[]
		for j,atom2 in enumerate(atoms):
			if atom2==main_elem and dists[i][j] < cutoff and (not j in done_atoms):
				vector = coords[j]-coords[i]
				#print(coords[i],coords2[i])
				unit = vector/np.linalg.norm(vector)
				coords_new = coords2[i]+unit*float(new_length)
				new_coords.append(list(coords_new))
				coords2[j]=coords_new
				new_atoms.append(main_elem)
				transfer_atoms.append(j)
				done_atoms.append(j)

ncoords = np.array(new_coords)
natoms = np.array(new_atoms)

new_ind_1 = natoms==main_elem
new_ind_2 = natoms==min_elem

new_dists, a, b, c, d = get_dists(ncoords, new_ind_1, new_ind_2)

new_elem_1_angles, new_elem_2_angles = get_angles(ncoords, new_ind_1, new_ind_2, new_dists, float(new_length)+0.01)

sorted_new=np.sort(new_elem_2_angles)
sorted_old =np.sort(elem_2_angles)

#print(sorted_new)
#print(sorted_old)

#good = []
#for i,atom in enumerate(sorted_old):
#	for j,angle in enumerate(atom):
#		good.append(math.isclose(angle,sorted_new[i][j]))

#print(good)

#if not (False in good):
#	print("It worked!")

write_xyz(name,natoms,ncoords)







