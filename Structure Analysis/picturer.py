import sys
from qchem_helper import *

xyz = sys.argv[1]

num_list = []
for num in sys.argv[2:]:
	num_list.append(int(num)+1) #+1 comes from -1 for 0 based indexing +2 for the first two lines of the .xyz

out_lines = []
with open(xyz, "r") as orig:
	for i, line in enumerate(orig):
		if i in num_list:
			out_lines.append(line)

with open("pics.xyz", "w") as write:
	write.write(""+str(len(out_lines))+" \n")
	write.write(" \n")
	for i, line in enumerate(out_lines):
		write.write(line)