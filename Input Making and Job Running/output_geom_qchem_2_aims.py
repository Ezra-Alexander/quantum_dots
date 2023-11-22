import numpy as np
import sys
from qchem_helper import get_geom_io

#from a qchem input or output (plot, not opt)

inpt = sys.argv[1]

atoms, coords = get_geom_io(inpt)

with open("geometry.in",'w') as aims:
	for i,atom in enumerate(atoms):
		aims.write("atom " + str(coords[i][0]) + " "+ str(coords[i][1]) + " "+ str(coords[i][2]) + " " + atom + "\n")