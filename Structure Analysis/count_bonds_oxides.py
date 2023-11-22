import numpy as np
import sys
from geom_helper import read_input_xyz, dist_all_points, dist_atom12, num_nn

def good_index(i, atoms, target):
	n_in = 0
	n_p = 0
	n_f = 0
	for j, atom in enumerate(atoms):
		if atom == "In" and target == "In":
			n_in=n_in+1
			if j==i:
				return n_in
		if atom == "P" and target == "P":
			n_p=n_p+1
			if j==i:
				return n_p
		if atom == "F" and target == "F":
			n_f=n_f+1
			if j==i:
				return n_f

def permutation(a,b):
    if len(a) != len(b):
        return False
    n = len(a)
    count = 0
    for i in range(n):
        for j in range(n):
            if b[j] == a[i]:
                count += 1
    return count == n


#to calculate the number of In-F bridging/terminal bonds from a qchem .xyz file

f_cutoff = 2.55
p_cutoff = 2.9
o_cutoff = 2.4
h_cutoff = 1
xyz = sys.argv[1]

coords, atoms = read_input_xyz(xyz)

index_f = atoms=="F"
index_p = atoms=="P"
index_in = atoms=="In"
index_o = atoms=="O"
index_h = atoms=="H"

dists = dist_all_points(coords)

#remember that this counts itself
f_nearest_neighbors = num_nn(dists, f_cutoff)
p_nearest_neighbors = num_nn(dists,p_cutoff)
o_nearest_neighbors = num_nn(dists, o_cutoff)
h_nearest_neighbors = num_nn(dists, h_cutoff)

n_terminal_f = 0
n_bridging_f = 0
n_triple_f=0

n_p_2 = 0
n_p_3 = 0
n_p_4 = 0

n_o_1 = 0
n_o_2 = 0
n_oh_1 = 0
n_oh_2 = 0


in_p_bonds = [-1]*len(atoms)

other = 0

#Code that helps determine appropriate cutoffs
for i, atom1 in enumerate(atoms):
	if index_in[i]:
		for j,atom2 in enumerate(atoms):
			if dists[i][j] > 2.55 and dists[i][j] < 3 and index_f[j]:
				print("we have a problem")
				print("distance", dists[i][j], "between F ", good_index(j,atoms,"F"), "and In ", good_index(i,atoms,"In"))
			if dists[i][j] > 2.9 and dists[i][j] < 3 and index_p[j]:
				print("we have a problem")
				print("distance", dists[i][j], "between P ", good_index(j,atoms,"P"), "and In ", good_index(i,atoms,"In"))
			if dists[i][j] > 2.4 and dists[i][j] < 3 and index_o[j]:
				print("we have a problem")
				print("distance", dists[i][j], "between O ", good_index(j,atoms,"O"), "and In ", good_index(i,atoms,"In"))

print()

for i, atom in enumerate(atoms):
	if index_f[i]:
		if f_nearest_neighbors[i] == 3:
			n_bridging_f += 1
		elif f_nearest_neighbors[i] == 2:
			n_terminal_f += 1
		elif f_nearest_neighbors[i] == 4:
			n_triple_f +=1
		else:
			print(good_index(i, atoms, "F"))
			other += 1
	if index_p[i]:
		if p_nearest_neighbors[i] == 3:
			n_p_2 += 1
		elif p_nearest_neighbors[i] == 4:
			n_p_3 += 1
		elif p_nearest_neighbors[i] == 5:
			n_p_4 += 1
		else:
			other += 1
	if index_o[i]:
		if o_nearest_neighbors[i]==4:
			n_oh_2+=1
		elif o_nearest_neighbors[i]==3:
			if np.count_nonzero(atoms=="H")>0:
				n_oh_1+=1
			else:
				n_o_2+=1
		else:
			other +=1

in_f_bonds = [-1]*len(atoms)
in_o_bonds = [-1]*len(atoms)

for i, atom in enumerate(atoms):
	if index_in[i]:
		in_f_bonds[i] = 0
		in_p_bonds[i] = 0
		in_o_bonds[i]=0
		for j,atom in enumerate(atoms):
			if dists[i][j] < 2.55 and index_f[j]:
				in_f_bonds[i]+=1
			if dists[i][j] < 2.9 and index_p[j]:
				in_p_bonds[i]+=1
			if dists[i][j] < 2.4 and index_o[j]:
				in_o_bonds[i]+=1


indiums = []

for i, atom in enumerate(atoms):
	if index_in[i]:
		indiums.append(str(in_p_bonds[i]+in_f_bonds[i]+in_o_bonds[i])+" , "+str(in_f_bonds[i])+ " , "+str(in_o_bonds[i]))

unique_in = set(indiums)

for i,ind in enumerate(unique_in):
	count = 0
	for j, atom in enumerate(indiums):
		if ind == indiums[j]:
			count += 1
	name = ind[:2] + "Coordinate Total" + ind[2:5] + " F" + ind[5:] + " Oxygen"
	print("There are", count, "Indium of type", name)
	print()

print("There are", n_p_2, "Two coordinate P,", n_p_3, "Three coordinate P, and", n_p_4, "Four coordinate P",)
print()

print("There are", n_triple_f, "Three-Coordinate F,", n_bridging_f, "bridging F and ", n_terminal_f, "terminal F in this structure.")
print()

print("There are", n_o_2, "bridging oxides, ", n_oh_1, "terminal hydroxides, and ", n_oh_2, "bridging hydroxides.")

if other != 0:
	print("you suck eggs")

#Try to find number of 4 member rings
rings = []
for i,atom in enumerate(atoms):
	for j,atom2 in enumerate(atoms):
		if dists[i][j] < 3.5 and i!=j:
			for k,atom3 in enumerate(atoms):
				if dists[i][k] < 3.5 and dists[j][k] < 3.5 and i!=k and j!=k:
					for l,atom4 in enumerate(atoms):
						if dists[i][l]<3.5 and dists[j][l]<3.5 and dists[k][l]<3.5 and i!=l and j!=l and k!=l:
							new = [i, j, k, l]
							count = 0
							for x in rings:
								if permutation(new,x):
									count = count+1
							if count==0:
								rings.append(new)

converted_rings = []
for ring in rings:
	count = 0
	for index in ring:
		if atoms[index] =="F" or atoms[index] =="O":
			count = count+1
	if count <3:
		converted_rings.append([atoms[ring[0]], good_index(ring[0], atoms, atoms[ring[0]]), atoms[ring[1]], good_index(ring[1], atoms, atoms[ring[1]]), atoms[ring[2]], good_index(ring[2], atoms, atoms[ring[2]]), atoms[ring[3]], good_index(ring[3], atoms, atoms[ring[3]])])

print()
print(converted_rings)
print()
print(len(converted_rings), "Four Membered Rings")


#[atom, good_index(i, atoms, atom), atom2, good_index(j, atoms, atom2), atom3, good_index(k, atoms, atom3), atom4, good_index(l, atoms, atom4)]
