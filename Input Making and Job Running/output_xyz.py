import numpy as np
import sys
from qchem_helper import get_converged_geom, write_xyz

#write an .xyz file for a converged output

outpt = sys.argv[1] #qchem.out

xyz_name=sys.argv[2]

coords,atoms = get_converged_geom(outpt)

write_xyz(xyz_name, atoms, coords)
