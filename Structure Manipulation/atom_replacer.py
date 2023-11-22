import sys
import numpy as np

#replaces all atoms of the first type with the second type in a given .xyz file

file = sys.argv[1]

to_replace = sys.argv[2]

replacememt = sys.argv[3]

with open(file, "r") as xyz:
	lines = xyz.readlines()
	for i, line in enumerate(lines):
		if line.find(to_replace)!= -1:
			lines[i] = replacememt + line[2:]


new_file = "new_" + file

with open(new_file, "w") as out:
	for i, line in enumerate(lines):
		out.write(line)