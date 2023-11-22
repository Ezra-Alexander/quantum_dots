import numpy as np
import sys
import math
from geom_helper import *
from qchem_helper import read_xyz, write_xyz

#a replacement of the geom_adder script
#this adds one of a specified atomic type to any number of atoms in your QD
#it tries to do this in the best way possible, including addition to 4c atoms

#The orginal .xyz file that I am editing
xyz = sys.argv[1]

#The name of the .xyz file this will make
name = sys.argv[2]

#what atom I want to add
add = sys.argv[3]

#the species you are adding to
add_base=sys.argv[4]

#the following inputs are any number of atom-specific
input_targets=[int(x) for x in sys.argv[5:]]

bond_max=3.0 #an arbitrary parameter for determining the upper distance between two atoms for them to be considered bonded
bond_lengths={"In":{"F":2.02,"Cl":2.37},"Ga":{"F":1.79,"Cl":2.19},"P":{"O":1.53}} #hard coded (terminal) bond lengths. not everything will be supported. Don't have to be perfect, these are just pre-opt estimates

coords,atoms=read_xyz(xyz)
temp_coords=np.copy(coords)
temp_atoms=np.copy(atoms)

#convert atom_specific to overall
count=0
targets=[]
for i,atom in enumerate(atoms):
	if atom==add_base:
		count=count+1
		if count in input_targets:
			targets.append(i)


#do the thing
for i,target in enumerate(targets):
	bond_length=bond_lengths[temp_atoms[target]][add]
	dists=dist_all_points(temp_coords)
	connectivity=connectivity_finder(dists,bond_max)

	if len(connectivity[target])==1:
		anchor=connectivity[target][0]
		anchor_bonds=connectivity[anchor]
		vector=temp_coords[anchor]-temp_coords[anchor_bonds[0]]
		unit=vector/np.linalg.norm(vector)
		added_coords=temp_coords[target]+(unit*bond_length)
		temp_coords=np.append(temp_coords,[added_coords],axis=0)
		temp_atoms=np.append(temp_atoms,add)

	elif len(connectivity[target])==2:
		v1=temp_coords[target]-temp_coords[connectivity[target][0]]
		v2=temp_coords[target]-temp_coords[connectivity[target][1]]

		parallel=(v1+v2)/np.linalg.norm(v1+v2)
		perp=np.cross(v1,v2)
		perp_unit=perp/np.linalg.norm(perp)

		adjacent=bond_length*math.cos(math.radians(109.5/2))
		opposite=bond_length*math.sin(math.radians(109.5/2))

		added_coords=temp_coords[target]+(parallel*adjacent)+(perp_unit*opposite)
		temp_coords=np.append(temp_coords,[added_coords],axis=0)
		temp_atoms=np.append(temp_atoms,add)

	elif len(connectivity[target])==3:
			temp_atoms,temp_coords=geom_adder(temp_atoms,temp_coords,target,add)

	elif len(connectivity[target])==4:
			#I think the best way to do this is to split the biggest gap, actually
			#but it has to be angle based
			max_angle=0 #arbitrary large number
			for j,bond1 in enumerate(connectivity[target]):
				for k,bond2 in enumerate(connectivity[target]):
					if j!=k:
						dist1 = dists[bond1][bond2]
						dist2 = dists[bond1][target]
						dist3 = dists[bond2][target]
						num = dist1**2 + dist2**2 - dist3**2
						denom = 2*dist1*dist2
						angle = math.acos(num/denom)
						if angle>max_angle:
							max_angle=angle
							max_index_1=bond1
							max_index_2=bond2

			v1=temp_coords[max_index_1]-temp_coords[target]
			v2=temp_coords[max_index_2]-temp_coords[target]

			vector=v1+v2
			unit=vector/np.linalg.norm(vector)

			added_coords=temp_coords[target]+(unit*bond_length)
			temp_coords=np.append(temp_coords,[added_coords],axis=0)
			temp_atoms=np.append(temp_atoms,add)



	else:
			print("you have a species here that is either 0 coordinate or more than 4 coordinate, neither of which is implemented yet")



write_xyz(name, temp_atoms, temp_coords)