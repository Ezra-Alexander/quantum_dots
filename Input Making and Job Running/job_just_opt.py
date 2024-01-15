import numpy as np
import sys
from qchem_helper import write_input, write_geom, my_get_rem_sp, read_xyz

#the .xyz file. Note that the file name is assumed to be unopt_DESCRIPTOR.xyz
xyz=sys.argv[1]
#template qchem opt .in
opt=sys.argv[2]
#charge
charge=sys.argv[3]

job_name=xyz[6:-4]

coords, atoms = read_xyz(xyz)

#the opt
opt_rem, trash = my_get_rem_sp(opt)

opt_rem =''.join(opt_rem)

write_input("opt_"+job_name+".in", atoms, coords, opt_rem, charge + " 1 \n")



