import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *

#a quick script that takes a .xyz and a specified atom (could I code this to take both index formats?)
#and prints the following geometry parameters: bond angles, bond lengths, coplanarity (3c only)

#formatted to be used to fill the "Defect Catalogue" excel spreadsheet

xyz=sys.argv[1] #the xyz file in question

atom_specific=False
if len(sys.argv) > 3: #for atom-specific indexing, i.e P 11
	target_atom = sys.argv[2]
	target_index = int(sys.argv[3])
	atom_specific = True
else: # for overall indexing, i.e. atom 111 
	overall_index = int(sys.argv[2])-1

cutoff = 2.9 #cutoff for 2 atoms to be considered bonded. Arbitrary and may need to be changed

coords,atoms = read_input_xyz(xyz)

dists=dist_all_points(coords)

#get all angles
angles=[]
for i,center in enumerate(coords):
    atom = []
    for j,atom1 in enumerate(coords):
        if dists[i][j]<cutoff and i!=j:
            for k,atom2 in enumerate(coords):
                if dists[i][k]<cutoff and i!=k and j!=k and k>j:
                    dist1 = dists[i][j]
                    dist2 = dists[i][k]
                    dist3 = dists[j][k]
                    num = dist1**2 + dist2**2 - dist3**2
                    denom = 2*dist1*dist2
                    angle = math.acos(num/denom)
                    atom.append(math.degrees(angle))
    angles.append(atom)

if atom_specific:
	target = 0
	for i,atom in enumerate(atoms):
		if atom == target_atom:
			target = target+1
			if target == target_index:
				overall_index = i

target_dists = [x for x in dists[overall_index] if 0<x<cutoff]
target_angles = angles[overall_index]
print()
print("Bond lengths are", end=" ")
[print(round(x,2), end=" ") for x in target_dists]
print()
print("Bond Angles are", end=" ")
[print(round(x,1), end=" ") for x in target_angles]
print()
print()

#let's compute, for 3c atoms only, the shortest distance between them and the plane made by their 3 bonding partners
#first we find the plane defined by the 3 bonding pairs

bonded_indeces = []
for i,atom in enumerate(atoms):
	if dists[overall_index][i] < cutoff and i != overall_index: # and atom!="F": #for computing In/Ga-4c
		bonded_indeces.append(i)

#bonded_indeces.pop(0) # for computing P-4c

if len(bonded_indeces) > 3:
	print("Hope this is supposed to be 4-coordinate:")

	min_coplane = 10 #arbitrary high number
	for i,bond in enumerate(bonded_indeces): #loop through bonded atoms, find coplane metric excluding that atom
		v1 = coords[bonded_indeces[i-1]] - coords[bonded_indeces[i-2]]
		v2 = coords[bonded_indeces[i-3]] - coords[bonded_indeces[i-2]]

		norm_v = np.cross(v1,v2) #plane equation comes from normal vector to plane

		a = norm_v[0]
		b = norm_v[1]
		c = norm_v[2]

		d = -a*coords[bonded_indeces[i-1]][0]-b*coords[bonded_indeces[i-1]][1]-c*coords[bonded_indeces[i-1]][2]

		dist_from_plane = (a*coords[overall_index][0]+b*coords[overall_index][1]+c*coords[overall_index][2]+d)/math.sqrt(a**2 + b**2 + c**2)

		if dist_from_plane<0:
			dist_from_plane = dist_from_plane*-1

		if dist_from_plane < min_coplane:
			min_coplane=dist_from_plane

	print("Most coplanar set of 3 has a distance from plane of",min_coplane)
	print()

elif len(bonded_indeces) < 3:
	print("Hope this is supposed to be 2-coordinate:")

	v1 = coords[bonded_indeces[1]] - coords[bonded_indeces[0]]
	u1 = v1/np.linalg.norm(v1)

	d1=coords[overall_index]-coords[bonded_indeces[0]]

	dist_vect=d1-np.dot(d1,u1)*u1

	dist_from_line=np.linalg.norm(dist_vect)

	print("Distance from center atom to line between 2 bonded atoms:",round(dist_from_line,2))
	print()

else:
	v1 = coords[bonded_indeces[1]] - coords[bonded_indeces[0]]
	v2 = coords[bonded_indeces[2]] - coords[bonded_indeces[0]]

	norm_v = np.cross(v1,v2) #plane equation comes from normal vector to plane

	a = norm_v[0]
	b = norm_v[1]
	c = norm_v[2]

	d = -a*coords[bonded_indeces[1]][0]-b*coords[bonded_indeces[1]][1]-c*coords[bonded_indeces[1]][2]

	dist_from_plane = (a*coords[overall_index][0]+b*coords[overall_index][1]+c*coords[overall_index][2]+d)/math.sqrt(a**2 + b**2 + c**2)

	if dist_from_plane<0:
		dist_from_plane = dist_from_plane*-1

	print("Distance from center atom to plane of 3 bonded atoms:",round(dist_from_plane,2))
	print()



