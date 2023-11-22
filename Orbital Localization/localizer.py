import numpy as np
import sys
from qchem_helper import write_input, write_geom, get_geom_io, my_get_rem_sp, get_geom_e_opt_last

#once again, the sbatch script is unfortunately hard-coded. 4 GB / core, 4 cores

#get ya inputs
input_file=sys.argv[1]
bas_file=sys.argv[2]
scr=sys.argv[3]
occ=sys.argv[4]
mo_indeces=sys.argv[5:]
mo_indeces=[int(x) for x in mo_indeces]

#the name of the localize input. Configured to work with "plot_job.in" type syntax
name = "loc_"+str(min(mo_indeces))+"_"+str(max(mo_indeces))+"_"+input_file[5:]

#get the index of the HOMO from nbas
nbas,nocc=np.loadtxt(bas_file,dtype=int,unpack=True)

#convert the mo_indeces to their more proper form
if occ=="o":
	for i,mo in enumerate(mo_indeces):
		mo_indeces[i]=nocc-mo
elif occ=="u":
	for i,mo in enumerate(mo_indeces):
		mo_indeces[i]=nocc+mo+1
else:
	print("Third variable should either be 'o' for occupied orbitals or 'u' for unoccupied orbitals")
	sys.exit()
mo_indeces.sort()

#open and read input file, inserting certain lines at the right moments
with open(input_file,"r") as inp:
	lines=[]
	rem_flag=False
	for i,line in enumerate(inp):
		if line.find("$rem") != -1: #edit 1: add rem variables for localization and auxilliary basis
			rem_flag=True
			lines.append(line)
		elif rem_flag:
			lines.append("ercalc 01104 \n")
			lines.append("aux_basis rij-def2-svp \n")
			lines.append("skip_scfman true \n")
			lines.append("scf_guess read \n")
			rem_flag=False
		elif line.find("alpha_molecular_orbital") != -1: #edit 2: only plot the orbitals we localize. Remember, localized orbitals are shifted to be occupied edge
			line = line[:24]+ str(nocc-len(mo_indeces)+1)+"-"+str(nocc) + " \n"
			lines.append(line)
		else:
			lines.append(line)

#edit 3: add $localize section
lines.append("\n")
lines.append("$localize \n")
for i, mo in enumerate(mo_indeces):
	lines.append(str(mo))
	lines.append(" ")
lines.append("\n")
lines.append("$end \n")

#write .in
with open(name,"w") as out:
	for i,line in enumerate(lines):
		out.write(line)


#write submit script
with open("submit_loc.sh",'w') as submit:
	submit.write("#!/bin/bash \n")
	submit.write("\n")
	submit.write("#SBATCH -J "+ name + " \n")
	submit.write("#SBATCH -o " + name[:-3] + ".log \n")
	submit.write("#SBATCH -e "+ name[:-3] + ".log \n")
	submit.write("#SBATCH --time unlimited \n")
	submit.write("#SBATCH -c 4 \n")
	submit.write("#SBATCH --mem-per-cpu 4000 \n")
	submit.write("\n")
	submit.write("scratch="+scr+" \n")
	submit.write("rm -r /scratch/ezraa/$scratch/ \n")
	submit.write("\n")
	submit.write("scp -r ../../$scratch/ . \n")
	submit.write("mv $scratch \"old_scratch\" \n")
	submit.write("\n")
	submit.write("scp -r \"old_scratch\" /scratch/ezraa \n")
	submit.write("\n")
	submit.write("qchem.latest -save -nt 4 " + name + " " + name[:-3] + ".out \"old_scratch\" \n")
	submit.write("\n")
	submit.write("cp -r $QCSCRATCH/\"old_scratch\" \"new_scratch\" \n")
	submit.write("\n")
	submit.write("rm -r /scratch/ezraa/\"old_scratch\" \n")
	submit.write("rm -r \"old_scratch\" \n")