import numpy as np
import sys
import math
from qchem_helper import read_xyz, write_xyz, write_input, get_rem_sp
from geom_helper import read_input_xyz, dist_all_points, dist_atom12, num_nn


#The orginal .xyz file that I am editing
xyz = sys.argv[1]

#The name of the .xyz file this will make
name = sys.argv[2]

coords, atoms = read_input_xyz(xyz)
dists = dist_all_points(coords)

with open(xyz,'r') as file:
	rem, spcharge = get_rem_sp(file.readlines())

indiums =[]
#loop over F atoms
#determine what Indium(s) each F is bonded to
for i, atom in enumerate(atoms):
	if atom =="F":
		bonds = []
		for j, atom2 in enumerate(atoms):
			if atom2=="In" and dists[i][j]< 2.6 :
				bonds.append(j)
		indiums.append(bonds)
	else:
		indiums.append(0)

#separate fluorines into bridging or terminal
bridging=[]
terminal=[]
for i,fluorine in enumerate(indiums):
	if type(fluorine) is list:
		if len(fluorine) ==2:
			bridging.append(True)
			terminal.append(False)
		if len(fluorine) ==1:
			bridging.append(False)
			terminal.append(True)
	else:
		bridging.append(False)
		terminal.append(False)

#Do terminal fluorines 
for i, fluorine in enumerate(indiums):
	if terminal[i]:
		indium = indiums[i][0]

		#make unit vector
		vector = coords[i]- coords[indium]
		mag = np.linalg.norm(vector)
		unit = vector/mag

		#replace F with O (or something else)
		oxygen_coords = coords[indium] + unit*2.017
		coords[i] = oxygen_coords
		atoms[i] = "Cl"

		#Add H to O
		#oh_dist_parallel = math.cos(math.radians(180-132))*1.145
		#normal_vector = np.cross(unit,[1, 0, 0]) / np.linalg.norm(np.cross(unit,[1, 1, 0]))
		#oh_dist_perp = math.sin(math.radians(180-132))*1.145
		#h_coords = oxygen_coords + unit*oh_dist_parallel + normal_vector*oh_dist_perp
		#coords=np.append(coords,[h_coords], axis=0)
		#atoms=np.append(atoms,"H")


#Do bridging fluorines
for i, fluorine in enumerate(indiums):
	if bridging[i]:
		indium1 = indiums[i][0]
		indium2 = indiums[i][1]

		#make unit vector
		vector = coords[i]- (coords[indium1]+coords[indium2])/2
		mag = np.linalg.norm(vector)
		unit = vector/mag

		#replace F with O
		oxygen_coords = (coords[indium1]+coords[indium2])/2 + unit*2.017
		coords[i] = oxygen_coords
		atoms[i] = "Cl"

		#Add H to O
		#oh_dist_parallel = math.cos(math.radians(180-132))*1.145
		#normal_vector = np.cross(unit,[1, 0, 0]) / np.linalg.norm(np.cross(unit,[1, 1, 0]))
		#oh_dist_perp = math.sin(math.radians(180-132))*1.145
		#h_coords = oxygen_coords + unit*oh_dist_parallel + normal_vector*oh_dist_perp
		#coords=np.append(coords,[h_coords], axis=0)
		#atoms=np.append(atoms,"H")

write_input(name,atoms,coords, rem, spcharge)
#write_xyz(name, atoms, coords)