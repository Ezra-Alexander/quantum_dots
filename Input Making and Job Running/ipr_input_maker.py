import numpy as np
import sys

#this should take a qchem plot input file and make a new double input file that gets what I need to do an IPR

input_file = sys.argv[1]
nbas = sys.argv[2]

with open(nbas,'r') as nbas:
	for i,line in enumerate(nbas):
		stripped = line.strip().split()
		basis = stripped[0]
		beta_electrons = stripped[1]
		n_virt = str(int(basis) - int(beta_electrons))

output_name = "ipr_" + input_file

original = []

new = []
new.append('\n')
new.append('@@@ \n')
new.append('\n')

charge_flag = False
rem_flag = False
plots_flag = False
with open(input_file,'r') as inp:
	for i,line in enumerate(inp):
		if line.find('$molecule') != -1:
			charge_flag=True
			new.append(line)
			original.append(line)
		elif charge_flag:
			c_m = line.strip().split()
			charge = c_m[0]
			mult = c_m[1]
			new_charge = str(int(charge) - 2*int(n_virt))	
			new.append(new_charge + " " + mult + "\n")
			charge_flag = False
			original.append(line)
		elif line.find('$rem') != -1:
			rem_flag = True
			new.append(line)
			original.append(line)
		elif rem_flag:
			new.append(line)
			new.append("lowdin_population 2 \n")
			new.append("scf_guess read \n")
			new.append("skip_scfman true \n")
			rem_flag = False
			original.append(line)
		elif line.find('plots true') != -1 or line.find('make_cube') != -1:
			print("No plots!")
		elif line.find('$plots') != -1:
			plots_flag = True
		elif line.find('$end') != -1 and plots_flag:
			plots_flag = False
		elif plots_flag:
			print("No plots!")
		else:
			new.append(line)
			original.append(line)

with open(output_name,'w') as file:
	for i,line in enumerate(original):
		file.write(line)

	for i, line in enumerate(new):
		file.write(line)




