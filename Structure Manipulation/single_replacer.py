import numpy as np
import sys
import math
from qchem_helper import read_xyz, write_xyz

#This only does very simple replacements - turn a single atom into an atom of a different type in the exact same location

#The orginal .xyz file that I am editing
xyz_in = sys.argv[1]

#The name of the .xyz file this will make
xyz_out = sys.argv[2]

#The index of the atom to change (total index, not atom-specific index) - converted to an int and made base 0 for python
target = int(sys.argv[3])-1

#what atom you want to turn it into
new_atom = sys.argv[4]

coords, atoms = read_xyz(xyz_in)

atoms[target]=new_atom

write_xyz(xyz_out, atoms, coords)