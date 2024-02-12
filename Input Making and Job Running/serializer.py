import sys
import numpy as np
from qchem_helper import *

#a quick script that takes a multi-.xyz file
#and writes a serialized qchem input file
#where each geometry is run as a single point (no plots)
#and the MOs from the previous geometry are used as the guess for the next geometry
#basis set can be read from a separate file

xyz=sys.argv[1] #the multi-.xyz. N Atoms \n comment \n atoms coords xN \n N Atoms ...
qchem=sys.argv[2] #name of the qchem file to write
basis = sys.argv[3] # either a .txt file containing a custom basis (no $basis) or the name of the basis to use in qchem format
functional = sys.argv[4] # the dft functional to use. qchem format
if len(sys.argv)>5:
	charge=sys.argv[5] #optional inputs for charge and spin. Otherwise assumes neutral singlet
	spin=sys.argv[6]
else:
	charge='0'
	spin='1'

with open(xyz,"r") as xyz_list:
	stucture_counter=-1
	structures=[]
	for i,line in enumerate(xyz_list):
		if i==0:
			n_atoms=int(line.strip().split()[0])
			xyz_length=n_atoms+2
		if i%xyz_length==0:
			stucture_counter=stucture_counter+1
			ith_structure=[]
		ith_structure.append(line)
		if (i+1)%xyz_length==0:
			structures.append(ith_structure)

n_files=len(structures)

basis_file=False
if basis[-4:]==".txt":
	basis_file=True

if basis_file:
	bas = open(basis, "r")
	bas_lines=bas.readlines()


with open(qchem,"w") as out:
	for i,structure in enumerate(structures):
		out.write("$molecule \n")
		out.write(charge+" "+spin+" \n")
		for j,line in enumerate(structure):
			if j>1:
				out.write(line)
		out.write('$end \n')
		out.write('\n')
		out.write("$rem \n")
		out.write('jobtype sp \n')
		out.write('method '+functional+' \n')
		if basis_file:
			out.write('basis GEN \n')
			out.write('purecart 1111 \n')
		else:
			out.write('basis '+basis+' \n')
		out.write('ithrsh_dft 15 \n')
		out.write('chelpg true\n')
		if i!=0:
			out.write('scf_guess read\n')
		out.write('$end \n')
		out.write('\n')
		if basis_file:
			out.write('$basis \n')
			for j,line in enumerate(bas_lines):
				out.write(line)
			out.write('$end \n')
			out.write('\n')
		if i!=len(structures)-1:
			out.write('@@@ \n')







