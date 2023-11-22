import numpy as np
import sys

#a script that adds the text of one file to the end of another file. designed for adding basis sets to an aims input, but it will work either way

inpt = sys.argv[1]

basis = sys.argv[2]

lines=[]
with open(inpt,'r') as inp:
	for i,line in enumerate(inp):
		lines.append(line)
	lines.append("\n")

with open(basis, 'r') as bs:
	for i, line in enumerate(bs):
		lines.append(line)

with open(inpt, 'w') as new:
	for i, line in enumerate(lines):
		new.write(line)