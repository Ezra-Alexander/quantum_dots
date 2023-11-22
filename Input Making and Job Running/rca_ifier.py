import numpy as np
import sys
from qchem_helper import write_input, write_geom, get_geom_io, my_get_rem_sp

#Takes a qchem input file of an opt that crashed and tries to make it better (algorithm to rca, increase max cycles, decrease threshold)

name = sys.argv[1]

atoms, coords = get_geom_io(name)

rem, sp = my_get_rem_sp(name)

rem.insert(1,"scf_algorithm rca \n")
rem.insert(1,"scf_convergence 6 \n")
rem.insert(1,"max_scf_cycles 500 \n")

new_rem =''.join(rem)

new_name = "2nd_" + name

write_input(new_name, atoms, coords, new_rem, sp)
