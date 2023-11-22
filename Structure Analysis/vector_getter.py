import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *

#returns a specified unit vector relative to some number of atoms in a structure

#parallel points towards the 1st input index

xyz=sys.argv[1] #.xyz file

vect=sys.argv[2] #specifies the vector you want. Either one perpendicular to the specified atoms ("perp") or one parallel ("para")

indices=[int(x)-1 for x in sys.argv[3:]] #overall indices please

coords,atoms = read_input_xyz(xyz)

output="Something went wrong"
if vect=="perp":
	if len(indices)!=3:
		print("3 points define a plane")
	else:
		v1=coords[indices[0]]-coords[indices[1]]
		v2=coords[indices[2]]-coords[indices[1]]

		normal=np.cross(v1,v2)
		output=normal/np.linalg.norm(normal)
elif vect=="para":
	if len(indices)!=2:
		print("2 points define a line")
	else:
		v1=coords[indices[0]]-coords[indices[1]]

		output=v1/np.linalg.norm(v1)
else:
	print("2nd argument should either be 'perp' for a vector perpendicular to those specified or 'para' for one parallel")

print()
print(output)
print()