import numpy as np
import sys
from geom_helper import write_xyz, get_geom_cp2k

inpt = sys.argv[1] #the cp2k multi-.xyz file. usually of the form CP2K-pos-1.xyz
out_name = sys.argv[2] #the name of the .xyz file to write

coords,atoms = get_geom_cp2k(inpt)

write_xyz(out_name, atoms, coords)