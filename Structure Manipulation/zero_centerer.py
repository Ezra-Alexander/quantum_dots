import numpy as np
import sys
import math
from qchem_helper import *

#this script takes a .xyz file or qchem input file of a structure not centered on (0,0,0) and generates the analagous structure centered on the origin 
#the center of mass of the new .xyz will be on the origin

xyz = sys.argv[1] #original file

name = sys.argv[2] #of the new file

#compute center of mass
#info the code has
masses = {"Ga":69.723, "P":30.974, "In":114.82,"Cl":35.453,"F":18.998403162}

if xyz[-4:]==".xyz":
	coords,atoms=read_xyz(xyz)
elif  xyz[-3:]==".in":
	atoms,coords=get_geom_io(xyz)
	rem,sp=my_get_rem_sp(xyz)
	new_rem=""
	for i,line in enumerate(rem):
		new_rem=new_rem+line
else:
	print("File type not supported")


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

#write the new file
if xyz[-4:]==".xyz":
	write_xyz(name, atoms, new_coords, comment='Shifted to be centered on the origin')
elif  xyz[-3:]==".in":
	write_input(name,atoms,new_coords,new_rem,sp,comment='Shifted to be centered on the origin')