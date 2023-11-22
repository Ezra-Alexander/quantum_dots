import numpy as np
import sys
import math
from geom_helper import *

#the sister script to defect_ifier.sh, meant to be called from that bash script but works otherwise

#takes an .xyz file and a number of indeces and makes a new .xyz file with the atoms with those indeces fully removed

#the .xyz file in question
xyz=sys.argv[1]

#the list of indeces of atoms to remove
indeces=sys.argv[2:]
indeces=[int(x) for x in indeces] 
indeces.sort(reverse=True)


#the part that matters
coords,atoms=read_input_xyz(xyz)
popped_atoms=[]
for i in indeces:
	coords=np.delete(coords,i-1,0)
	popped_atoms.append(atoms[i-1])
	atoms=np.delete(atoms,i-1,0)


#coming up with a name
atom_string=""
for atom in popped_atoms:
	atom_string=atom_string+atom+"_"

atom_string=atom_string[:-1]

name="unopt_"+atom_string+xyz[5:]


write_xyz(name,atoms,coords)