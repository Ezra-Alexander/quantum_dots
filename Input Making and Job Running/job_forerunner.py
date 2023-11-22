import numpy as np
import sys
from qchem_helper import write_input, write_geom, my_get_rem_sp, read_xyz

#the .xyz file, without the ".xyz"
xyz=sys.argv[1]
#template qchem opt .in
#opt=sys.argv[2]
#template qchem single point .in
sp=sys.argv[2]
#charge
charge=sys.argv[3]


coords, atoms = read_xyz(xyz+".xyz")


#the opt
#opt_rem, trash = my_get_rem_sp(opt)

#opt_rem =''.join(opt_rem)

#write_input("opt_"+xyz+".in", atoms, coords, opt_rem, charge + " 1 \n")



#the single point

sp_rem, trash = my_get_rem_sp(sp)

#count electrons - right now this only works for the elements In, P, F, O, and H (but it's easy to add more)

n_elec = 0
for i,atom in enumerate(atoms):
	if atom=="H":
		n_elec = n_elec+1
	if atom=="O":
		n_elec = n_elec+8
	if atom=="F":
		n_elec = n_elec+9
	if atom=="P":
		n_elec=n_elec+15
	if atom=="In":
		n_elec=n_elec+21
	if atom=="Li":
		n_elec = n_elec+3

n_elec = n_elec//2

sp_rem.append('\n')
sp_rem.append( '$plots \n')
sp_rem.append( 'grid_range (-20,20) (-20,20) (-20,20) \n')
sp_rem.append( 'grid_points 150 150 150 \n')
sp_rem.append( "alpha_molecular_orbital " + str(n_elec-9) + "-" + str(n_elec+10) + " \n")
sp_rem.append( "$end \n")

sp_rem =''.join(sp_rem)

write_input("plot_"+xyz+".in", atoms, coords, sp_rem, charge + " 1 \n")



#submit script

with open("submit_plot.sh",'w') as submit:
	submit.write("#!/bin/bash \n")
	submit.write("\n")
	submit.write("#SBATCH -J "+ "plot_"+xyz+".in" + " \n")
	submit.write("#SBATCH -o " + "plot_"+xyz+".log" + "\n")
	submit.write("#SBATCH -e "+ "plot_"+xyz+".log" + "\n")
	submit.write("#SBATCH --time unlimited \n")
	submit.write("#SBATCH -c 4 \n")
	#submit.write("#SBATCH --mem-per-cpu 4000 \n")
	submit.write("\n")
	submit.write("scratch=\"scratch_"+ "plot_"+xyz + "\" \n")
	submit.write("curr_d=$PWD \n")
	submit.write("\n")
	submit.write("rm -r $QCSCRATCH/$scratch \n")
	submit.write("cp -r $scratch $QCSCRATCH \n")
	submit.write("\n")
	submit.write("qchem.latest -save -nt 4 " + "plot_"+xyz+".in" + " " + "plot_"+xyz+".out" + " $scratch \n")
	submit.write("\n")
	submit.write("cp -r $QCSCRATCH/$scratch $curr_d \n")




