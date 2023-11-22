import numpy as np
import sys
from qchem_helper import get_geom_io

#.xyz only

inpt = sys.argv[1]

with open(inpt,'r') as inp:
	geom = []
	for i,line in enumerate(inp):
		if i > 1:
			geom.append(line.strip().split())

	geom = np.array(geom)
	atoms = geom[:,0]
	coords = geom[:,1:].astype(float)

with open("geometry_2.in",'w') as aims:
	for i,atom in enumerate(atoms):
		aims.write("atom " + str(coords[i][0]) + " "+ str(coords[i][1]) + " "+ str(coords[i][2]) + " " + atom + "\n")