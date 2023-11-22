import numpy as np
import sys
import math
from geom_helper import *

#for cutting up the output of interpolations

xyz=sys.argv[1]

with open(xyz,"r") as inp:
	length=0
	total=0
	lines=[]
	for i,line in enumerate(inp):
		if i==0:
			length=int(line.strip().split()[0])+2
		total=total+1
		lines.append(line)

for i in range(total//length):
	with open(str(i)+".xyz","w") as out:
		for j,line in enumerate(lines):
			if j>=(i*length) and j<((i+1)*length):
				out.write(line)
