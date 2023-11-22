import numpy as np
import sys
from qchem_helper import write_input, write_geom, get_geom_io, my_get_rem_sp, get_geom_e_opt_last

#Takes a qchem input file of an opt that crashed and tries to make it better (algorithm to rca, increase max cycles, decrease threshold)

inpt = sys.argv[1]
oupt = sys.argv[2]

out_lines = open(oupt, "r")
geom = get_geom_e_opt_last(out_lines.readlines())

atoms = geom[:,0]
coords = geom[:,1:].astype(float)

rem, sp = my_get_rem_sp(inpt)

rem.insert(1,"geom_opt_update 5 \n")
rem.insert(1,"geom_opt_max_cycles 500 \n")


new_rem =''.join(rem)

new_name = "restart_" + inpt

write_input(new_name, atoms, coords, new_rem, sp)

