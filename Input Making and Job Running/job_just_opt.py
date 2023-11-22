import numpy as np
import sys
from qchem_helper import write_input, write_geom, my_get_rem_sp, read_xyz

#the .xyz file, without the ".xyz"
xyz=sys.argv[1]
#template qchem opt .in
opt=sys.argv[2]
#charge
charge=sys.argv[3]


coords, atoms = read_xyz(xyz+".xyz")


#the opt
opt_rem, trash = my_get_rem_sp(opt)

opt_rem =''.join(opt_rem)

write_input("opt_"+xyz+".in", atoms, coords, opt_rem, charge + " 1 \n")



