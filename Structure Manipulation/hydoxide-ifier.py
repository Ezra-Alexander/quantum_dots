import numpy as np
import sys
import math
from qchem_helper import read_xyz, write_xyz

#I'm not sure this will work for bridging hydroxides but whatever

#The orginal .xyz file that I am editing
xyz = sys.argv[1]

#The name of the .xyz file this will make
name = sys.argv[2]

#The index of the indium the F is on (total index, not In index) - converted to an int and made base 0 for python
indium = int(sys.argv[3])-1
#The index of the F to be turned into Hydroxide (total index, not F index) - converted to an int and made base 0 for python
fluorine = int(sys.argv[4])-1

coords, atoms = read_xyz(xyz)

#make unit vector
vector = coords[fluorine]- coords[indium]
mag = np.linalg.norm(vector)
unit = vector/mag

#replace F with O
oxygen_coords = coords[indium] + unit*2.017
coords[fluorine] = oxygen_coords
atoms[fluorine] = "O"

#Add H to O
oh_dist_parallel = math.cos(math.radians(180-132))*1.145
normal_vector = np.cross(unit,[1, 0, 0]) / np.linalg.norm(np.cross(unit,[1, 0, 0]))
oh_dist_perp = math.sin(math.radians(180-132))*1.145
h_coords = oxygen_coords + unit*oh_dist_parallel + normal_vector*oh_dist_perp
coords=np.append(coords,[h_coords], axis=0)
atoms=np.append(atoms,"H")


write_xyz(name, atoms, coords)