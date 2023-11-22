import numpy as np
import sys
import math
from geom_helper import *

# a script for growing InLx layers on top of a P-3c terminated facet

xyz=sys.argv[1] #the .xyz file
to_write=sys.argv[2] #name of the new .xyz file to write
ligand=sys.argv[3] #ligand species
growth_type=sys.argv[4] #to use a 3 P-3c facet as an example, you can grow either In3Cl6 or In3Cl6P, depending on the orientation of the In
#use P to indicate the P-centered type, and noP to indicate the other type

#the following inputs should be at least one P-specific index for each facet you want to grow on
facets_in=[[int(x)] for x in sys.argv[5:]]

bond_max=2.83 #a manual parameter for bond lengths


#let's start by determining the complete facets
coords,atoms=read_input_xyz(xyz)
dists=dist_all_points(coords)
connectivity=connectivity_finder(dists,bond_max)

#convert to overall indices, adjusted for python
overall_facets=[[] for x in facets_in]
cations=[]
for i,facet in enumerate(facets_in):
	count=0
	for j,atom in enumerate(atoms):
		if atom=="P":
			count=count+1
			if count==facet[0]:
				overall_facets[i]=[j]
		if atom=="In" or atom=="Ga" or atom=="Al":
			cations.append(atom)

#complete the facets
for i,facet in enumerate(overall_facets):
	for j,cation in enumerate(connectivity[facet[0]]):
		for k,bond in enumerate(connectivity[cation]):
			if len(connectivity[bond])==3 and atoms[bond]=="P" and (bond not in facet):
				overall_facets[i].append(bond)

#add an In (or Ga) to each P
temp_atoms=np.copy(atoms)
temp_coords=np.copy(coords)
new_in=[[] for x in overall_facets]
for i,facet in enumerate(overall_facets):
	for j,phos in enumerate(facet):
		#lift cation to add from .xyz
		new_atoms,new_coords=geom_adder(temp_atoms,temp_coords,phos,cations[((i)*len(facet))+j])
		temp_atoms=np.copy(new_atoms)
		temp_coords=np.copy(new_coords)
		new_in[i].append(len(temp_atoms)-1)


#Fill in the new layer, depending on the value of 'growth_type'

if growth_type=="P": #THIS IS NOT DONE YET
	for i,facet in enumerate(new_in):
		p_set=[] #the set of indium coordinates. i don't know why i called it that
		for j,indium in enumerate(facet): #grab the In on each facet
			p_set.append(new_coords[indium])

		if len(p_set)==3:
			p_set=np.array(p_set)
			p_flip=np.transpose(p_set)
			p_xy=np.average(p_flip,axis=1)

			v1=p_set[0]-p_set[1] #which way is up?
			v2=p_set[0]-p_set[2]
			up=np.cross(v1,v2)
			up=up/np.linalg.norm(up)

			ref_bond=p_set[0]-new_coords[overall_facets[i][0]] #gotta correct for the "randomness" of the cross product
			direction=math.degrees(math.acos(np.dot(up,ref_bond)/(np.linalg.norm(ref_bond)*np.linalg.norm(up))))
			if direction>90:
				up=-1*up

			angle=120
			p_coords=np.copy(p_xy)
			while angle>109.5: #this is such a bad way of doing this
				p_coords=p_coords+(up*0.1)
				v1=new_coords[new_in[i][0]]-p_coords
				v2=new_coords[new_in[i][1]]-p_coords
				angle=math.degrees(math.acos((np.dot(v1,v2))/(np.linalg.norm(v1)*np.linalg.norm(v2))))
			new_atoms=np.append(new_atoms,"P")
			new_coords=np.append(new_coords,[p_coords],axis=0)

			#Now I need to add 2 F to each In/Ga

			for j,indium in enumerate(facet):
				new_atoms, new_coords = bring_to_4c(new_atoms,new_coords,indium,ligand,3.5)

			

		elif len(p_set)==1:
			new_atoms, new_coords = bring_to_4c(new_atoms,new_coords,facet[0],ligand,bond_max)


		elif len(p_set)>3:
			print("P set facet growth not supported yet for facets larger than 3 P")
		else:
			print("Need at least 3 P-3c to grow on")


elif growth_type=="noP":

	for i,facet in enumerate(new_in):
		p_set=[]
		for j,indium in enumerate(facet): #grab the In on each facet
			p_set.append(new_coords[indium])

		if len(p_set)==1:
			new_atoms, new_coords = bring_to_4c(new_atoms,new_coords,facet[0],ligand,bond_max)

		else:

			print("F set facet growth not coded yet")

else:
	print("Accepted inputs for the 4th argument are 'P' and 'noP'. Check script notes for more details")


write_xyz(to_write, new_atoms, new_coords)


