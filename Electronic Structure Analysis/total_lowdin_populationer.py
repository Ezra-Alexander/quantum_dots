import numpy as np
import sys
import math

#the goal is to write a script that can go through a variety of qchem output files where Lowdin populations have been printed and sum them up for various atoms
#in an ipr.out type double file, would like to be able to do both ground state and specified unoccupied states

#for now let's write this to work for just a single file, I can mass script-ify it later

out_file=sys.argv[1]
#n_virt=sys.argv[2] #if the file is of the ipr double file format, this is the number of virtual orbitals to include

#hard coding number of protons, accounting for effective core potentials
n_protons={"In":21,"Ga":31,"Cl":17,"F":9,"P":15,"Al":13}


#start reading the file
with open(out_file,"r") as out:
	job=0
	total_lowdin=[]
	total_lowdin_2=[]
	low_read=0

	for i,line in enumerate(out):
		if  line.find("day")!=-1:
			job=job+1
		elif line.find("beta electrons")!=-1 and job==0:
			n_occ=int(line.strip().split()[2])
		elif line.find("beta electrons")!=-1 and job==1:
			n_occ_2=int(line.strip().split()[2])
		elif line.find("Partial Lowdin")!=-1:
			low_read=1
		elif low_read==1:
			low_read=2
		elif low_read>1 and line.find("Atomic Becke Charges")!=-1:
			low_read=0
		elif low_read>1 and job==0:
			total_lowdin.append(line.strip().split())
		elif low_read>1 and job==1:
			total_lowdin_2.append(line.strip().split())

#get the lowdin arrays into a nice format
total_lowdin.pop(-1)
total_lowdin.pop(-1)

#counting
total_length=len(total_lowdin)
count=0
species_list=[]
for i,line in enumerate(total_lowdin):
	line_length=len(line)
	if i==0:
		short_length=line_length
	else:
		if line_length==short_length:
			n_rows=count
			break
		else:
			count=count+1
			species_list.append(line[1])

#pop off the indexing lines
npops=math.floor(total_length/n_rows)
for i in range(npops):
	backwards=npops-i-1
	total_lowdin.pop((n_rows+1)*backwards)

#line them up
orderly_lowdin = [[] for i in range(n_rows)]
remainder=n_occ%6
if remainder==0:
	remainder=6
new_length=len(total_lowdin)
for i,line in enumerate(total_lowdin):
	target_row=int(line[0])-1
	if i>=(new_length-n_rows):
		orderly_lowdin[target_row].extend(line[(-1*remainder):])
	else:
		orderly_lowdin[target_row].extend(line[-6:])

orderly_lowdin=np.array(orderly_lowdin)
orderly_lowdin=orderly_lowdin.astype(float)

#print
orbs=["s","p","d"]
if(n_rows%3)!=0:
	print("You have a problem: code assumed exactly s,p,d orbitals for all atoms")
for i,species in enumerate(species_list):
	print(species,"     ",orbs[i%3],"     ", np.sum(orderly_lowdin[i]))
	# if (i%3)==2:
	# 	print(species,"     ","total","     ", -(np.sum(orderly_lowdin[i])+np.sum(orderly_lowdin[i-1])+np.sum(orderly_lowdin[i-2]))+n_protons[species[:2]])



#if 2 jobs, do it again
if job==2:
	total_lowdin_2.pop(-1)
	total_lowdin_2.pop(-1)

	total_length_2=len(total_lowdin_2)


