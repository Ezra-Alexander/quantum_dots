import sys
from qchem_helper import get_geom_io

#note that rows (nshells) doesn't work if there is H in the structure

output_file=sys.argv[1]

atoms, coords = get_geom_io(output_file)
nshells = 3*len(atoms)

basis = 0
with open(output_file,'r') as out:
	for i,line in enumerate(out):
		if line.find('basis functions') != -1:
			basis = line.strip().split()[5]
			break

with open("rows_basis.txt", "w") as txt:
	txt.write(str(nshells)+" "+basis)