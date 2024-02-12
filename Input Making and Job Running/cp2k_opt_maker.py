import numpy as np
import sys
from geom_helper import *

#given a .xyz and reference cp2k.inp and run_cp2k.sh
#makes a cp2k.inp with the right geometry
#makes a run_cp2k.sh with the right file name
#doesn't allow for the finagling of parameters in either file (reads everything from references)

xyz=sys.argv[1] #.xyz file
ref_opt=sys.argv[2] #reference cp2k.inp opt
ref_sub=sys.argv[3] #reference submit script

#read_everything
coords,atoms=read_input_xyz(xyz)

opt=open(ref_opt,"r")
opt_lines=opt.readlines()

sub=open(ref_sub,"r")
sub_lines=sub.readlines()

flag=0
descriptor=""
for i,char in enumerate(xyz):
	if char=="_" and flag==0:
		flag=1
	elif char==".":
		flag=0
	elif flag==1:
		descriptor=descriptor+char

#write opt
with open("opt_cp2k.inp","w") as inp:
	geom_flag=0
	for i,line in enumerate(opt_lines):
		if line.find("&COORD") != -1:
			geom_flag=geom_flag+1
			inp.write(line)
		elif geom_flag>0 and line.find("&END COORD") != -1:
			if geom_flag<=len(atoms): #handles case where reference geometry section is too short
				n_more=len(atoms)-geom_flag+1
				for j in range(n_more):
					inp.write(atoms[geom_flag-1+j]+" "+str(coords[geom_flag-1+j][0])+" "+str(coords[geom_flag-1+j][1])+" "+str(coords[geom_flag-1+j][2])+" \n")
			geom_flag=0
			inp.write(line)
		elif geom_flag>0:
			if geom_flag<=len(atoms):
				inp.write(atoms[geom_flag-1]+" "+str(coords[geom_flag-1][0])+" "+str(coords[geom_flag-1][1])+" "+str(coords[geom_flag-1][2])+" \n")
			geom_flag=geom_flag+1
		else:
			inp.write(line)

#write submit
with open("run_cp2k.sh","w") as sub:
	for i,line in enumerate(sub_lines):
		if line.find("--job-name")!= -1:
			sub.write("#SBATCH --job-name=opt_"+descriptor+" \n")
		elif line.find("--output")!= -1:
			sub.write("#SBATCH --output=opt_cp2k.out \n")

		elif line.find('srun')!= -1:
			srun_params=line.strip().split()
			flag=False
			for j,term in enumerate(srun_params):
				if flag:
					sub.write("opt_cp2k.inp")
					sub.write(" \n")
				else:
					sub.write(term)
					sub.write(" ")
					if term=="-i":
						flag=True

		else:
			sub.write(line)

