import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *

#a variation on the geom_param_getter.py script that repeats for several atoms and then prints averages
#only configured to work with atom specific


xyz=sys.argv[1] #the xyz file in question
element=sys.argv[2] #the element of all specified atoms. should be the same
indices=[int(x) for x in sys.argv[3:]] #just a list of element specific indices (don't repeat the elements)

cutoff = 3.0 #cutoff for 2 atoms to be considered bonded. Arbitrary and may need to be changed

#info the code has
masses = {"Ga":69.723, "P":30.974, "In":114.82,"Cl":35.453,"F":18.998403162,"Al":26.982, "Zn":65.38,"S":32.06, "Se":78.971,"Si":28.085}

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

target = 0
overall_indices=[] #convert easy to input indices into overall indices
for i,atom in enumerate(atoms):
	if atom == element:
		target = target+1
		if target in indices:
			overall_indices.append(i)

max_list=[]
for i, atom in enumerate(dists):
	max_list.append(max(atom))
print(max(max_list))


#choose center of dot
#compute center of mass
cent =np.array([0.0])
mass = 0
for i,atom in enumerate(atoms):
	cent = cent + coords[i]*masses[atom]
	mass = mass + masses[atom]
	com = cent/mass

#do the shift
new_coords = coords-com


dist_set=[] #get the distances and angles and coplanarity and zero-centered coordinates for all of our targets
ang_set=[]
coplane_set=[]
coords_set=[]
for i,ind in enumerate(overall_indices):
	#print(get_atom_specific_index(atoms,ind))
	target_dists = [x for x in dists[ind] if 0<x<cutoff]
	#print(target_dists)
	#print(target_dists)
	target_dists.sort()
	target_angles = angles[ind]
	target_angles.sort()
	#if max(target_dists)>2.75:
	#	print(get_atom_specific_index(atoms,ind))
	#if max(target_angles)>133:
	#	print(get_atom_specific_index(atoms,ind))
	#if min(target_angles)<189 and max(target_angles)>129:
	#	print(get_atom_specific_index(atoms,ind))
	coplanarity=get_coplane(coords,dists,ind,cutoff)
	#print(coplanarity)
	#if coplanarity < 0.45:
	#	print(get_atom_specific_index(atoms,ind))
	target_coords = new_coords[ind]

	dist_set.append(target_dists)
	ang_set.append(target_angles)
	coplane_set.append(coplanarity)
	coords_set.append(target_coords)

dist_set=np.array(dist_set)
ang_set=np.array(ang_set)
coplane_set=np.array(coplane_set)
coords_set=np.array(coords_set)

dist_means=np.mean(dist_set,axis=0)
ang_means=np.mean(ang_set,axis=0)
coplane_mean=np.mean(coplane_set)
coords_means=np.mean(coords_set,axis=0)
print()
print("Mean bond lengths are", end=" ")
for i in dist_means:
	print(round(i,2), end=" ")
print()
print()
print("Mean bond angles are", end=" ")
for i in ang_means:
	print(round(i,1), end=" ")
print()
print()
print("Mean coplanarity metric is",round(coplane_mean,2))
print()
print("Mean coordinate vector (COM adjusted) is", end=" ")
for i in coords_means:
	print(round(i,2), end=" ")
print()
print()

print("max bond length is:", np.max(dist_set, axis=0))
print("min bond length is:", np.min(dist_set, axis=0))
print("max bond angle is:", np.max(ang_set, axis=0))
print("min bond angle is:", np.min(ang_set, axis=0))
print("average bond angle is:", np.average(ang_set))


