import numpy as np
import sys
from qchem_helper import *

#the sister script to the bash script of the same name
#meant to take in a .xyz file and some parameters and turn them into a qchem single point input file and submit script
#the plan is to try combining the ipr and plot jobs, but i'm not sure if that will actually work with the scratch directory yet

xyz=sys.argv[1]#the .xyz file. Hopefully the only file input
cores=sys.argv[2]#the number of cores you want the calculation to run on
plot_min=int(sys.argv[3]) #the lowest energy MO you want to plot. 0 for HOMO, 1 for HOMO-1, etc
plot_max=int(sys.argv[4]) #the highest energy MO you want to plot. 0 for LUMO, 1 for LUMO+1, etc
mode=sys.argv[5] #just a plot ("plot") or the ipr+plot ("both")
mem_static=sys.argv[6] #the mem_static value in qchem's $rem, in MB
priority=sys.argv[7] #the priority for the job. If on Ulysses, say 'ulysses'. priotity automatically sets mem_total

telemachus=True
if priority=='ulysses':
	telemachus=False

#disctionaries with built-in info
atomic_charges= {"In":3, "Ga":3, "P":-3, "F":-1, "Cl":-1, "Al":3,"Zn":2,"Se":-1,"S":-2,"O":-2,"H":1} #manual dictionary of the preferred oxidation state of each atom
n_basis_functions = {'In': 26, 'Ga': 32, 'P': 18, 'Cl': 18, 'Br': 32,'F':14,'H':5,'O':14, 'C':14,'S':18, 'Li':9, 'Al':18, "Zn":31,"S":18,"Se":32, "Si":18} #assuming def2-svp
n_electrons = {'In': 21, 'Ga': 31, 'P': 15, 'Cl': 17, 'Br': 35,'F':9,'H':1,'O':8, 'C':14,'S':6, 'Li':3, 'Al':13, "Zn":30,"S":16,"Se":34, "Si":14} #assuming def2-ecp. For neutral atom

#read .xyz
coords, atoms = read_xyz(xyz)

# re-center coordinates around (0,0,0)
average=[0, 0, 0]
for i,atom in enumerate(coords):
	average=average+atom
average=average/len(atoms)
coords=coords-average

#get electronic descriptors
charge=0
nbas=0
nelec=0
for i,atom in enumerate(atoms):
	charge=charge+atomic_charges[atom]
	nbas=nbas+n_basis_functions[atom]
	nelec=nelec+n_electrons[atom]-atomic_charges[atom]
nocc=int(nelec/2)
spch_string=str(charge)+" 1"

#determine MO bounds for plot
lower_bound=nocc-plot_min
upper_bound=nocc+1+plot_max
bounds=str(lower_bound)+"-"+str(upper_bound)

#if doing the joint plot+ipr
if mode=="both":

	#make file name
	descriptor=xyz[6:-4]
	file_name="spipr_"+descriptor+".in"


	#determine ipr "charge"
	n_virt = nbas - nocc
	new_charge = str(charge - 2*n_virt)

	#write .in
	if telemachus:
		write_ipr(file_name,atoms,coords,spch_string,new_charge,static=mem_static,bound_string=bounds,cluster=priority)
	else:
		write_ipr(file_name,atoms,coords,spch_string,new_charge,static=mem_static,bound_string=bounds)

elif mode=="plot":

	#make file name
	descriptor=xyz[6:-4]
	file_name="sp_"+descriptor+".in"

	#write .in
	if telemachus:
		my_write_plot(file_name,atoms,coords,spch_string,static=mem_static,bound_string=bounds,cluster=priority)
	else:
		my_write_plot(file_name,atoms,coords,spch_string,static=mem_static,bound_string=bounds)

elif mode=="ipr":

	#make file name
	descriptor=xyz[6:-4]
	file_name="ipr_"+descriptor+".in"

	#determine ipr "charge"
	n_virt = nbas - nocc
	new_charge = str(charge - 2*n_virt)

	#write .in
	if telemachus:
		write_ipr(file_name,atoms,coords,spch_string,new_charge,static=mem_static,cluster=priority)
	else:
		write_ipr(file_name,atoms,coords,spch_string,new_charge,static=mem_static)

else:
	raise Exception("Chosen mode not yet implemented!")


#write submit_plot.sh
if telemachus:
	write_submit(file_name,cores,mode,cluster=priority)
else:
	write_submit(file_name,cores,mode)	

