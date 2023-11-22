import numpy as np
import sys
import math
from geom_helper import *

#this code is meant to take all bonds of a given type (i.e P-In bonds) in a given structure and adjust them to a new bond length
#as written this only works for two element systems
#doesn't work with 4-membered rings

#Different from the original in two ways
#works (exclusively) on double .xyzs as input and output
#actually changes the majority element to Lithium

#i want to edit this such that relative changes in bond length from the ideal are preserved
#note that for now this only works for In-P

xyz = sys.argv[1] #original .xyz file

name = sys.argv[2] #of the new xyz file

element1 = sys.argv[3] #the test is coded such that the minority element is element 2 for some reason
element2 = sys.argv[4]

new_length = float(sys.argv[5]) #the ideal bond length

cutoff = 2.9 #set manually, may need to be changed. #make sure it is greater than the minimum distance between identical elements

ideal_inp=2.579 #set manually, from the bulk crystal

with open(xyz,"r") as inp:
	length=0
	atoms_start=[]
	coords_start=[]
	atoms_end=[]
	coords_end=[]
	for i,line in enumerate(inp):
		if i==0:
			length=int(line.strip().split()[0])+2
		elif i==1:
			print("Converting "+element1+" to Li")
		elif i < length:
			entries = line.strip().split()
			atoms_start.append(entries[0])
			coords_start.append([float(entries[1]),float(entries[2]), float(entries[3])])
		elif i>length+1:
			entries = line.strip().split()
			atoms_end.append(entries[0])
			coords_end.append([float(entries[1]),float(entries[2]), float(entries[3])])


atoms_start=np.array(atoms_start)
coords_start=np.array(coords_start)
atoms_end=np.array(atoms_end)
coords_end=np.array(coords_end)


#do the start
ind_1 = atoms_start==element1
ind_2 = atoms_start==element2

dists, a, b, c, d = get_dists(coords_start, ind_1, ind_2)

elem_1_angles, elem_2_angles = get_angles(coords_start, ind_1, ind_2, dists, cutoff)

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
coords2 = coords_start.copy()

for i, atom in enumerate(atoms_start):
	if atom == min_elem and count==0:
		new_coords[0]=list(coords_start[i])
		new_atoms.append(min_elem)
		count = count+1
		done_atoms.append(i)

		for j,atom2 in enumerate(atoms_start):
			if atom2==main_elem and dists[i][j] < cutoff:
				vector = coords_start[j]-coords_start[i]
				old_length=np.linalg.norm(vector)
				unit = vector/old_length
				coords_new = coords_start[i]+(unit*new_length*(old_length/ideal_inp))
				new_coords.append(list(coords_new))
				coords2[j]=coords_new
				new_atoms.append(main_elem)
				transfer_atoms.append(j)
				done_atoms.append(j)

	elif atom == min_elem and count>0:
		
		count=count+1
		for j,atom2 in enumerate(atoms_start):
			if atom2==main_elem and dists[i][j] < cutoff and j in transfer_atoms:				
				vector = coords_start[i]-coords_start[j]
				old_length=np.linalg.norm(vector)
				unit = vector/old_length
				coords_new = coords2[j]+(unit*new_length*(old_length/ideal_inp))
				new_coords.append(list(coords_new))
				coords2[i]=coords_new
				new_atoms.append(min_elem)
				done_atoms.append(i)

		transfer_atoms=[]
		for j,atom2 in enumerate(atoms_start):
			if atom2==main_elem and dists[i][j] < cutoff and (not j in done_atoms):
				vector = coords_start[j]-coords_start[i]
				old_length=np.linalg.norm(vector)
				#print(coords[i],coords2[i])
				unit = vector/old_length
				coords_new = coords2[i]+(unit*new_length*(old_length/ideal_inp))
				new_coords.append(list(coords_new))
				coords2[j]=coords_new
				new_atoms.append(main_elem)
				transfer_atoms.append(j)
				done_atoms.append(j)

ncoords_start = np.array(new_coords)
natoms_start = np.array(new_atoms)


#now repeat for the second .xyz

ind_1 = atoms_end==element1
ind_2 = atoms_end==element2

dists, a, b, c, d = get_dists(coords_end, ind_1, ind_2)

elem_1_angles, elem_2_angles = get_angles(coords_end, ind_1, ind_2, dists, cutoff)

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
coords2 = coords_end.copy()

for i, atom in enumerate(atoms_end):
	if atom == min_elem and count==0:
		new_coords[0]=list(coords_end[i])
		new_atoms.append(min_elem)
		count = count+1
		done_atoms.append(i)

		for j,atom2 in enumerate(atoms_end):
			if atom2==main_elem and dists[i][j] < cutoff:
				vector = coords_end[j]-coords_end[i]
				old_length=np.linalg.norm(vector)
				unit = vector/old_length
				coords_new = coords_end[i]+(unit*new_length*(old_length/ideal_inp))
				new_coords.append(list(coords_new))
				coords2[j]=coords_new
				new_atoms.append(main_elem)
				transfer_atoms.append(j)
				done_atoms.append(j)

	elif atom == min_elem and count>0:
		
		count=count+1
		for j,atom2 in enumerate(atoms_end):
			if atom2==main_elem and dists[i][j] < cutoff and j in transfer_atoms:				
				vector = coords_end[i]-coords_end[j]
				old_length=np.linalg.norm(vector)
				unit = vector/old_length
				coords_new = coords2[j]+(unit*new_length*(old_length/ideal_inp))
				new_coords.append(list(coords_new))
				coords2[i]=coords_new
				new_atoms.append(min_elem)
				done_atoms.append(i)

		transfer_atoms=[]
		for j,atom2 in enumerate(atoms_end):
			if atom2==main_elem and dists[i][j] < cutoff and (not j in done_atoms):
				vector = coords_end[j]-coords_end[i]
				old_length=np.linalg.norm(vector)
				#print(coords[i],coords2[i])
				unit = vector/old_length
				coords_new = coords2[i]+(unit*new_length*(old_length/ideal_inp))
				new_coords.append(list(coords_new))
				coords2[j]=coords_new
				new_atoms.append(main_elem)
				transfer_atoms.append(j)
				done_atoms.append(j)

ncoords_end = np.array(new_coords)
natoms_end = np.array(new_atoms)


#some tests to check if it worked (the tests don't work, which is why the're commented out. Just check manually)
#new_ind_1 = natoms==main_elem
#new_ind_2 = natoms==min_elem
#new_dists, a, b, c, d = get_dists(ncoords, new_ind_1, new_ind_2)
#new_elem_1_angles, new_elem_2_angles = get_angles(ncoords, new_ind_1, new_ind_2, new_dists, float(new_length)+0.01)
#sorted_new=np.sort(new_elem_2_angles)
#sorted_old =np.sort(elem_2_angles)
#print(sorted_new)
#print(sorted_old)
#good = []
#for i,atom in enumerate(sorted_old):
#	for j,angle in enumerate(atom):
#		good.append(math.isclose(angle,sorted_new[i][j]))
#print(good)
#if not (False in good):
#	print("It worked!")


#now we need to change all the cations to Lithium
li_atoms_start=[]
for i,atom in enumerate(natoms_start):
	if atom==element1:
		li_atoms_start.append("Li")
	else:
		li_atoms_start.append(atom)

li_atoms_end=[]
for i,atom in enumerate(natoms_end):
	if atom==element1:
		li_atoms_end.append("Li")
	else:
		li_atoms_end.append(atom)

#now we need to stitch the two back together

with open(name,'w') as xyz_file:
	xyz_file.write(str(len(li_atoms_start))+'\n')
	xyz_file.write("Pristine \n")
	for i,atom in enumerate(li_atoms_start):
		xyz_file.write(atom)
		xyz_file.write("  ")
		xyz_file.write(str(ncoords_start[i][0]))
		xyz_file.write("  ")
		xyz_file.write(str(ncoords_start[i][1]))
		xyz_file.write("  ")
		xyz_file.write(str(ncoords_start[i][2]))
		xyz_file.write(" \n")
	xyz_file.write(str(len(li_atoms_end))+'\n')
	xyz_file.write("Pristine \n")
	for i,atom in enumerate(li_atoms_end):
		xyz_file.write(atom)
		xyz_file.write("  ")
		xyz_file.write(str(ncoords_end[i][0]))
		xyz_file.write("  ")
		xyz_file.write(str(ncoords_end[i][1]))
		xyz_file.write("  ")
		xyz_file.write(str(ncoords_end[i][2]))
		xyz_file.write(" \n")






