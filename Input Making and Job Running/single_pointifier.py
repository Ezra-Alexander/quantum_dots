import numpy as np
import sys
from qchem_helper import write_input, write_geom, get_geom_io, my_get_rem_sp, get_geom_e_opt_last

#inpt is the name of the optimization input file 
inpt = sys.argv[1]
oupt = sys.argv[2]

with open(oupt,'r') as out:
        flag = 0
        flag2 = 1
        geom = []
        for i,line in enumerate(out):
        	if line.find('beta') != -1 and flag2 == 1:
        		nelec = int(line.strip().split()[2])
        		flag2 = 0
        	if line.find('Z-matrix Print') != -1 and flag >= 1:
        		break
        	if flag > 1:
        		geom.append(line.strip().split()[1:])
        	if flag > 0:
        		flag += 1
        	if line.find('OPTIMIZATION CONVERGED') != -1:
        		flag = 1

geom = np.array(geom[3:-1])
atoms = geom[:,0]
coords = geom[:,1:].astype(float)


rem, sp = my_get_rem_sp(inpt)


#this is hard coded
new_rem = ['$rem\n', 'jobtype sp  \n', 'method pbe0 \n', 'basis def2-svp \n', 'ecp def2-ecp \n', 'mem_total 4000 \n', 'mem_static 1000 \n', 'ithrsh_dft 15 \n', 'make_cube_files true \n', 'plots true\n', '$end \n', '\n', '$plots \n', 'grid_range (-20,20) (-20,20) (-20,20) \n', 'grid_points 150 150 150 \n']

nelec_line = "alpha_molecular_orbital " + str(nelec-9) + "-" + str(nelec+10) + " \n"

new_rem.append(nelec_line)
new_rem.append("$end \n")
new_rem =''.join(new_rem)


test = True

while test:
	if inpt[0] == "2":
		inpt = inpt[4:]
	elif inpt[0] == "o":
		inpt = inpt[4:]
	else:
		test = False

new_name = "plot_optzd_" + inpt

write_input(new_name, atoms, coords, new_rem, sp)


with open("new_submit.sh",'w') as submit:
	submit.write("#!/bin/bash \n")
	submit.write("\n")
	submit.write("#SBATCH -J "+ new_name + " \n")
	submit.write("#SBATCH -o " + new_name[:-3] + ".log \n")
	submit.write("#SBATCH -e "+ new_name[:-3] + ".log \n")
	submit.write("#SBATCH --time unlimited \n")
	submit.write("#SBATCH -c 4 \n")
	submit.write("#SBATCH --mem-per-cpu 4000 \n")
	submit.write("\n")
	submit.write("scratch=\"scratch_"+ new_name[:-3] + "\" \n")
	submit.write("curr_d=$PWD \n")
	submit.write("\n")
	submit.write("rm -r $QCSCRATCH/$scratch \n")
	submit.write("cp -r $scratch $QCSCRATCH \n")
	submit.write("\n")
	submit.write("qchem.latest -save -nt 4 " + new_name + " " + new_name[:-3] + ".out $scratch \n")
	submit.write("\n")
	submit.write("cp -r $QCSCRATCH/$scratch $curr_d \n")









