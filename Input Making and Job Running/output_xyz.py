import numpy as np
import sys
from qchem_helper import get_converged_geom, write_xyz

#write an .xyz file for a converged output

outpt = sys.argv[1]

xyz_name=outpt[:-4]+".xyz"

atoms, coords = get_converged_geom(outpt)

write_xyz(xyz_name, atoms, coords)
