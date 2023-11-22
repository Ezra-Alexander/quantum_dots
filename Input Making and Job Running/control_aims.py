import numpy as np
import sys

#this assumes we aren't changing the charge, multiplicity, functional, basis set, etc
#all we really need to change is the range on compute_dipolematrix and output dos/species_proj_dos
#for compute_dipolematrix we want E_(HOMO+2) - 1 eV to E_(LUMO+2) + 1 eV and for the dos it's +/- 3
#which means we need to read in the qchem output

#still need to sort out what needs to be in the directory and what doesn't

ref=sys.argv[1]

qchem=sys.argv[2]


#step 1 - get the E_HOMO+2 and E_LUMO+2


with open(qchem,'r') as q:
	flag = 0
	mo_energies = []
	for i,line in enumerate(q):
		if line.find(' -- Virtual --') != -1 and flag == 1:
			break
		elif line.find('-- Occupied --') != -1 and flag == 0:
			flag=1
		elif flag ==1:
			mo_energies.append(line.strip().split())

occ_energies=[]
for energy in mo_energies:
	occ_energies += energy

with open(qchem,'r') as q:
	flag = 0
	mo_energies = []
	for i,line in enumerate(q):
		if line.find('--------------------------------------------------------------') != -1 and flag == 1:
			break
		elif line.find('-- Virtual --') != -1 and flag == 0:
			flag=1
		elif flag ==1:
			mo_energies.append(line.strip().split())

virt_energies = []
for energy in mo_energies:
	virt_energies += energy

#don't forget to convert to eV!
e_homo_2 = 27.2114*float(occ_energies[-3])
e_lumo_2 = 27.2114*float(virt_energies[2])


#step 2 - edit control.in to have the correct bounds 
#again, this does not currently change charge, multiplicity, anything


with open(ref,'r') as r:
	lines = []
	for i,line in enumerate(r):
		if line.find('compute_dipolematrix') != -1 and line.find('#')==-1:
			lines.append('compute_dipolematrix ' + str(e_homo_2 - 1) + " "+ str(e_lumo_2+1) +" 1 \n")
		elif line.find('output dos') != -1 and line.find('#')==-1:
			lines.append('output dos ' + str(e_homo_2 - 3) + " "+ str(e_lumo_2+3) +" 2000 0.1 \n")
		elif line.find('output species_proj_dos') != -1 and line.find('#')==-1:
			lines.append('output species_proj_dos ' + str(e_homo_2 - 3) + " "+ str(e_lumo_2+3)+ " 2000 0.1 \n")
		else:
			lines.append(line)

with open("control.in", "w") as cont:
	for i, line in enumerate(lines):
		cont.write(line)


#step 3 - edit submit_aims.sh so that the job name is different


with open("submit_aims.sh",'r') as sub:
	lines = []
	for i,line in enumerate(sub):
		if line.find(" -J ") != -1:
			lines.append("#SBATCH -J " + "aims_dipole"+ ref[5:-4] + " \n")
		else:
			lines.append(line)

with open("submit_aims.sh", "w") as submit:
	for i, line in enumerate(lines):
		submit.write(line)