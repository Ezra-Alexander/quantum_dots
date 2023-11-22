import numpy as np
import sys
from qchem_helper import write_xyz, get_geom_io

inpt = sys.argv[1]

atoms, coords = get_geom_io(inpt)

write_xyz("optzd_"+inpt[5:-3]+".xyz", atoms, coords)
