import numpy as np
import sys
import math
from geom_helper import *

pristine=sys.argv[1]
defect=sys.argv[2]
p_indeces=sys.argv[3:]
p_indeces=[int(x) for x in p_indeces]

threshold=3.0 #somewhat arbitrary, to determine if a P is bonded to an In. May need to be changed for certain systems

#find all In bound to indexed P in each dot

pris_coords, pris_atoms = read_input_xyz(pristine)
def_coords, def_atoms=read_input_xyz(defect)

#pristine first
#pull out overall indeces of target P
pris_overall=[]
p_count = 0
for i,atom in enumerate(pris_atoms):
	if atom=="P":
		p_count=p_count+1
		if p_count in p_indeces:
			pris_overall.append(i)

pris_In = pris_atoms=='In'
pris_P = pris_atoms=='P'
pris_Ga = pris_atoms=="Ga"
pris_InGa = np.logical_or(pris_In,pris_Ga)
pris_lig = np.logical_not(np.logical_or(pris_InGa,pris_P))

all_pris_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(pris_coords,pris_InGa,pris_P,pris_lig) #don't care about other dists, just an artefact of the helper I'm using

pris_in_overall=[]
for i,atom in enumerate(pris_atoms):
	if pris_InGa[i]:
		for j,phos in enumerate(pris_overall):
			if all_pris_dists[i][phos]<threshold:
				if i not in pris_in_overall:
					pris_in_overall.append(i)


#now repeat it all for the defected dot
#pull out overall indeces of target P
def_overall=[]
p_count = 0
for i,atom in enumerate(def_atoms):
	if atom=="P":
		p_count=p_count+1
		if p_count in p_indeces:
			def_overall.append(i)

def_In = def_atoms=='In'
def_P = def_atoms=='P'
def_Ga = def_atoms=="Ga"
def_InGa = np.logical_or(def_In,def_Ga)
def_lig = np.logical_not(np.logical_or(def_InGa,def_P))

all_def_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(def_coords,def_InGa,def_P,def_lig) #don't care about other dists, just an artefact of the helper I'm using

def_in_overall=[]
for i,atom in enumerate(def_atoms):
	if def_InGa[i]:
		for j,phos in enumerate(def_overall):
			if all_def_dists[i][phos]<threshold:
				if i not in def_in_overall:
					def_in_overall.append(i)


#check to make sure everything is ok
if len(def_in_overall)!= len(pris_in_overall):
	print("Problem! The bonding has changed in the defected dot, and the endpoints have different # of In")

n_atoms=len(def_overall)+len(def_in_overall)

#now, write endpoints.xyz
with open("endpoints.xyz","w") as out:
	out.write(str(n_atoms)+"\n")
	out.write("Pristine \n")
	for i,index in enumerate(pris_overall):
		out.write("P")
		out.write("  ")
		out.write(str(pris_coords[index][0]))
		out.write("  ")
		out.write(str(pris_coords[index][1]))
		out.write("  ")
		out.write(str(pris_coords[index][2]))
		out.write(" \n")
	for i,index in enumerate(pris_in_overall):
		out.write(str(pris_atoms[index]))
		out.write("  ")
		out.write(str(pris_coords[index][0]))
		out.write("  ")
		out.write(str(pris_coords[index][1]))
		out.write("  ")
		out.write(str(pris_coords[index][2]))
		out.write(" \n")
	out.write(str(n_atoms)+"\n")
	out.write("Defect \n")
	for i,index in enumerate(def_overall):
		out.write("P")
		out.write("  ")
		out.write(str(def_coords[index][0]))
		out.write("  ")
		out.write(str(def_coords[index][1]))
		out.write("  ")
		out.write(str(def_coords[index][2]))
		out.write(" \n")
	for i,index in enumerate(def_in_overall):
		out.write(str(def_atoms[index]))
		out.write("  ")
		out.write(str(def_coords[index][0]))
		out.write("  ")
		out.write(str(def_coords[index][1]))
		out.write("  ")
		out.write(str(def_coords[index][2]))
		out.write(" \n")
