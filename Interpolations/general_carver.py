import numpy as np
import sys
import math
from geom_helper import *

#i'm going to assume (for now) that, if you are using the bulk crystal as a reference, you are only targeting one atom
#i'm now going to work on removing that assumption

reference=sys.argv[1] #.xyz
final_structure=sys.argv[2] #.xyz
target_species=sys.argv[3] #
target_indeces=sys.argv[4:] #atom specific, any number
target_indeces=[int(x) for x in target_indeces]

threshold=2.9 #somewhat arbitrary, to determine if two atoms are bonded. May need to be changed for certain systems

#read xyzs
pris_coords, pris_atoms = read_input_xyz(reference)
def_coords, def_atoms=read_input_xyz(final_structure)

#for the final structure
#extract the targeted atoms
def_targets=[]
count = 0
for i,atom in enumerate(def_atoms):
	if atom==target_species:
		count=count+1
		if count in target_indeces:
			def_targets.append(i)

#get all bond distances
def_In = def_atoms=='In'
def_P = def_atoms=='P'
def_Ga = def_atoms=="Ga"
def_InGa = np.logical_or(def_In,def_Ga)
def_lig = np.logical_not(np.logical_or(def_InGa,def_P))
lig="F" #I give myself a chance to overwrite this later

all_def_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(def_coords,def_InGa,def_P,def_lig) #don't care about other dists, just an artefact of the helper I'm using

#now pull out everything bonded to your targets
def_bonded=[]
bonded_atom_order=[]
ligand_count=0 #how many ligands are bound to our target. for interpolating bulk crystal In to surface In w/ ligands. structured as a list, where each entry is the number of F on each target
for i,atom in enumerate(def_atoms):
	if atom!=target_species:
		for j,targ in enumerate(def_targets):
			if all_def_dists[i][targ]<threshold:
				if i not in def_bonded:
					def_bonded.append(i)
					bonded_atom_order.append(atom)
					if def_lig[i]:
						ligand_count=ligand_count+1
						lig=atom


#now, if your target is In/Ga, you need one extra layer of Li to charge balance
#adding H to P
layer_count=0
natoms_third_layer=0
if target_species=="In" or target_species=="Ga":
	third_layer=[]
	for j,bond in enumerate(def_bonded):
		if def_atoms[bond]=="P":
			sub_set=[]
			layer_count=layer_count+1
			for i,atom in enumerate(def_atoms):
				if all_def_dists[i][bond]<threshold:
					if (i not in def_bonded) and (i not in def_targets):

				#		#need to prevent double-counting
				#		double_counting=False
				#		if layer_count>1:
				#			for k in range(layer_count):
				#				if i in third_layer[k-1]:
				#					double_counting=True 
				#		if not double_counting:

							sub_set.append(i)
							natoms_third_layer=natoms_third_layer+1

			third_layer.append(sub_set)


#now for the intial structure

#get all bond distances
pris_In = pris_atoms=='In'
pris_P = pris_atoms=='P'
pris_Ga = pris_atoms=="Ga"
pris_InGa = np.logical_or(pris_In,pris_Ga)
pris_lig = np.logical_not(np.logical_or(pris_InGa,pris_P))

all_pris_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(pris_coords,pris_InGa,pris_P,pris_lig) #don't care about other dists, just an artefact of the helper I'm using

#extract the targeted atoms
pris_targets=[]
count = 0
if target_species=="In":
	target_indeces[0]=666 #this is shameless hardcoding for In, but it makes my life so much easier. Will need to find Ga and P analogues eventually
if target_species=="Ga":
	target_indeces[0]=378 #this is shameless hardcoding for In, but it makes my life so much easier. Will need to find Ga and P analogues eventually
if target_species=="P" and reference=="InP.xyz":
	target_indeces[0]=669

for i,atom in enumerate(pris_atoms):
	if atom==target_species:
		count=count+1
		if count in target_indeces:
			pris_targets.append(i)

#these are specific to both initial and final structures. don't @ me
if len(target_indeces)==2: #more blatant hard-coding. don't forget these are 0-based indices
	if target_species=="P": #assumes InP.xyz reference
		pris_targets[0]=1340
		pris_targets[1]=1269
	if target_species=="Ga":
		pris_targets[0]=377
		pris_targets[1]=2664

#now pull out everything bonded to your targets
pris_bonded=[]
pris_lig_count=0
for i,atom in enumerate(pris_atoms):
	if atom!=target_species:
		for j,targ in enumerate(pris_targets):
			if all_pris_dists[i][targ]<threshold:
				if i not in pris_bonded:
					pris_bonded.append(i)
					if pris_lig[i]:
						pris_lig_count=pris_lig_count+1


#one remaining problem
#the initial atoms all need to be as close as possible to their final counterparts
#num_flips=0
for j,final_atom in enumerate(def_bonded):
	current_dist=np.linalg.norm(pris_coords[pris_bonded[j]]-def_coords[final_atom]-pris_coords[pris_targets[0]]+def_coords[def_targets[0]])
	min_dist=1000
	min_i=10	
	for i,inital_atom in enumerate(pris_bonded):
		distance=np.linalg.norm(pris_coords[inital_atom]-def_coords[final_atom]-pris_coords[pris_targets[0]]+def_coords[def_targets[0]])
		if i!=j and distance<min_dist:
			min_dist=distance
			min_i=i
	if j!=min_i and min_dist<current_dist: #and num_flips<1:
#		print(j,def_atoms[final_atom], current_dist)
#		print(min_i,pris_atoms[pris_bonded[min_i]],min_dist)

		#temp_coords=pris_coords[pris_bonded[j]]
		#print(temp_coords)
		#temp_atoms=pris_atoms[pris_bonded[j]]

		#pris_coords[pris_bonded[j]]=pris_coords[pris_bonded[min_i]]
		#pris_atoms[pris_bonded[j]]=pris_atoms[pris_bonded[min_i]]

		#print(temp_coords)
		#pris_coords[pris_bonded[min_i]]=temp_coords
		#pris_atoms[pris_bonded[min_i]]=temp_atoms

		#num_flips=num_flips+1

		pris_coords[[pris_bonded[j],pris_bonded[min_i]]]=pris_coords[[pris_bonded[min_i], pris_bonded[j]]]


#now I may need to substitute some pristine bonded atoms to ligand atoms
if pris_lig_count<ligand_count:

		#some futzing around with bond lengths
		bond_length=1.986 #terminal In-F bond length
		ideal_inp_bond=2.579 #from bulk crystal
		if target_species=="Ga":
				bond_length=1.774 #terminal Ga-F bond length
		if lig!="F":
			print("Error - You need to hard code this cation-ligand bond length")

		new_order=[]
		for i,atom in enumerate(bonded_atom_order):
			if atom!="P":
				pris_atoms[pris_bonded[i]]=lig

				#if there is more than one target, i need to determine the closer target
				min_dist=1000
				closer_target_index=pris_targets[0]
				bridging_flag=False
				if len(target_indeces)>1:
					for j, index in enumerate(pris_targets):
						distance = np.linalg.norm(pris_coords[pris_bonded[i]]-pris_coords[index])
						if abs(distance-min_dist)<.2:
							bridging_flag=True
						if distance < min_dist:
							min_dist=distance
							closer_target_j=j

				bond=pris_coords[pris_bonded[i]]-pris_coords[pris_targets[closer_target_j]]
				reference_bond_length=np.linalg.norm(def_coords[def_bonded[i]]-def_coords[def_targets[closer_target_j]])
				normal=bond/np.linalg.norm(bond)
				new_coords=pris_coords[pris_targets[closer_target_j]]+(normal*reference_bond_length)
				#new_coords=pris_coords[pris_targets[closer_target_j]]+(normal*bond_length)
				#if bridging_flag and len(target_indeces)==2:
				#	other_bond=pris_coords[pris_bonded[i]]-pris_coords[pris_targets[closer_target_j-1]]
				#	other_normal=other_bond/np.linalg.norm(other_bond)
				#	other_coords=pris_coords[pris_targets[closer_target_j-1]]+(other_normal*reference_bond_length)
				#	#other_coords=pris_coords[pris_targets[closer_target_j-1]]+(other_normal*bond_length)
				#	new_coords= (new_coords+other_coords)/2
				pris_coords[pris_bonded[i]]=new_coords
			
elif pris_lig_count>ligand_count:
	print("boy what happened you have more initial ligand than final ligand")
	


#check to make sure our number of bonded atoms hasn't changed
if len(def_bonded)!= len(pris_bonded):
	print("Problem! The bonding has changed in the defected dot, and the endpoints have different # of atoms")
	print(len(def_bonded), len(pris_bonded))

n_atoms=len(def_targets)+len(def_bonded)

#gotta recompute the pristine distances after the swaps
all_pris_dists, inp_dists, pin_dists, inf_dists, inpf_dists = get_dists(pris_coords,pris_InGa,pris_P,pris_lig) #don't care about other dists, just an artefact of the helper I'm using


#now we need to implement the third layer in both dots
if target_species=="In" or target_species=="Ga":
	#we assume that there are no In-P-In-P 4-member rings here	
	#lip_bond_length=2.384
	lip_bond_length=1.421 #actually H-P bond length ...

	#first we do the hydrogen substitutions in the final structure
	count=0
	added=[]
	for i,phos in enumerate(def_bonded):
		if def_atoms[phos]=="P":
			for j, li in enumerate(third_layer[count]):
				#def_atoms[li]="Li"
				if li not in added:
					def_atoms[li]="H" #lithium dominates the Cb edge, trying protons
					bond=def_coords[li]-def_coords[phos]
					normal=bond/np.linalg.norm(bond)
					new_coords=def_coords[phos]+(lip_bond_length*normal)
					def_coords[li]=new_coords
					added.append(li)
				else:
					def_atoms=np.append(def_atoms,"H")
					bond=def_coords[li]-def_coords[phos]
					normal=bond/np.linalg.norm(bond)
					new_coords=def_coords[phos]+(lip_bond_length*normal)
					def_coords=np.append(def_coords,[new_coords], axis=0)
					third_layer[count][j]=len(def_atoms)-1
			count=count+1

	#when you have multiple target In/Ga, and bridging atoms in the 3rd layer, and double up on H, you have a steric clash
	#i'm going to see if I can resolve this by rotating half the PH3 groups
	to_rotate=[]
	for i,phos_1 in enumerate(third_layer):
		for j,h_1 in enumerate(phos_1):
			for k,phos_2 in enumerate(third_layer):
				if i<k: #prevents double-counting
					for l,h_2 in enumerate(phos_2):
						dist=np.linalg.norm(def_coords[h_2]-def_coords[h_1])
						if dist<1.75:
							to_rotate.append(k)
	

	#now I need to rotate all H on the target P by ~60 degrees
	original_def_coords=def_coords.copy()
	for i, phos in enumerate(to_rotate):
		#find the central P
		p_count=0
		for k, bonded in enumerate(def_bonded):
			if def_atoms[bonded]=="P":
				if p_count==phos:
					center_coords=def_coords[bonded]
					break
				else:
					p_count=p_count+1

		for j, proton in enumerate(third_layer[phos]):

			#this is almost certainly not the best way to do this
			proton_1_coords=def_coords[proton]
			proton_2_coords=original_def_coords[third_layer[phos][j-1]]
			average_coords=(proton_2_coords+proton_1_coords)/2
			bond=average_coords-center_coords
			normal=bond/np.linalg.norm(bond)
			new_coords=center_coords+(normal*lip_bond_length)
			def_coords[proton]=new_coords



	#now we compute the initial 3rd layer
	pris_third_layer=[]
	for j,bond in enumerate(pris_bonded):
		if pris_atoms[bond]=="P":
			sub_set=[]
			for i,atom in enumerate(pris_atoms):
				if all_pris_dists[i][bond]<threshold:
					if (i not in pris_bonded) and (i not in pris_targets):
						sub_set.append(i)
			pris_third_layer.append(sub_set)

	#now we must add H to the initial structure
	count=0
	added=[]
	for i,phos in enumerate(pris_bonded):
		if pris_atoms[phos]=="P":
			for j, li in enumerate(pris_third_layer[count]):
				#pris_atoms[li]="Li"
				if li not in added:
					pris_atoms[li]="H"
					bond=pris_coords[li]-pris_coords[phos]
					normal=bond/np.linalg.norm(bond)
					new_coords=pris_coords[phos]+(lip_bond_length*normal)
					pris_coords[li]=new_coords
					added.append(li)
				else:
					pris_atoms=np.append(pris_atoms,"H")
					bond=pris_coords[li]-pris_coords[phos]
					normal=bond/np.linalg.norm(bond)
					new_coords=pris_coords[phos]+(lip_bond_length*normal)
					pris_coords=np.append(pris_coords,[new_coords], axis=0)
					pris_third_layer[count][j]=len(pris_atoms)-1
			count=count+1

	#repeating for pristone
	#when you have multiple target In/Ga, and bridging atoms in the 3rd layer, and double up on H, you have a steric clash
	#i'm going to see if I can resolve this by rotating half the PH3 groups
	to_rotate=[]
	for i,phos_1 in enumerate(pris_third_layer):
		for j,h_1 in enumerate(phos_1):
			for k,phos_2 in enumerate(pris_third_layer):
				if i<k: #prevents double-counting
					for l,h_2 in enumerate(phos_2):
						dist=np.linalg.norm(pris_coords[h_2]-pris_coords[h_1])
						if dist<1.75:
							to_rotate.append(k)


	#now I need to rotate all H on the target P by ~60 degrees in pristine too
	original_pris_coords=pris_coords.copy()
	for i, phos in enumerate(to_rotate):
		#find the central P
		p_count=0
		for k, bonded in enumerate(pris_bonded):
			if pris_atoms[bonded]=="P":
				if p_count==phos:
					center_coords=pris_coords[bonded]
					break
				else:
					p_count=p_count+1
		for j, proton in enumerate(pris_third_layer[phos]):
				#this is almost certainly not the best way to do this
			proton_1_coords=pris_coords[proton]
			proton_2_coords=original_pris_coords[pris_third_layer[phos][j-1]]
			average_coords=(proton_2_coords+proton_1_coords)/2
			bond=average_coords-center_coords
			normal=bond/np.linalg.norm(bond)
			new_coords=center_coords+(normal*lip_bond_length)
			pris_coords[proton]=new_coords


	n_atoms=len(def_targets)+len(def_bonded)+natoms_third_layer


	#now we have to swap the ordering of the hydrogen in the initial structure to match the final structure
	for i,final_subset in enumerate(third_layer):
		for j, final_h_index in enumerate(final_subset):
			current_dist=np.linalg.norm(pris_coords[pris_third_layer[i][j]]-def_coords[final_h_index]-pris_coords[pris_targets[0]]+def_coords[def_targets[0]])
			min_dist=1000
			min_i=10	
			for k,inital_h_index in enumerate(pris_third_layer[i]):
				distance=np.linalg.norm(pris_coords[inital_h_index]-def_coords[final_h_index]-pris_coords[pris_targets[0]]+def_coords[def_targets[0]])
				if k!=j and distance<min_dist:
					min_dist=distance
					min_k=k
			if j!=min_i and min_dist<current_dist:
				pris_coords[[pris_third_layer[i][j],pris_third_layer[i][min_k]]]=pris_coords[[pris_third_layer[i][min_k], pris_third_layer[i][j]]]


#now, write endpoints.xyz
with open("endpoints.xyz","w") as out:
	out.write(str(n_atoms)+"\n")
	out.write("Pristine \n")

	for i,index in enumerate(pris_targets):
		if i==1: #these if elses are temporary cheese on my part
			out.write(target_species)
			out.write("  ")
			out.write(str(pris_coords[index][0]+.392))
			out.write("  ")
			out.write(str(pris_coords[index][1]))
			out.write("  ")
			out.write(str(pris_coords[index][2]-.392))
			out.write(" \n")
		else:
			out.write(target_species)
			out.write("  ")
			out.write(str(pris_coords[index][0]))
			out.write("  ")
			out.write(str(pris_coords[index][1]))
			out.write("  ")
			out.write(str(pris_coords[index][2]))
			out.write(" \n")

	count=0
	for i,index in enumerate(pris_bonded):
		if i==0 or i==3 or i==5:
			out.write(str(pris_atoms[index]))
			out.write("  ")
			out.write(str(pris_coords[index][0]+.392))
			out.write("  ")
			out.write(str(pris_coords[index][1]))
			out.write("  ")
			out.write(str(pris_coords[index][2]-.392))
			out.write(" \n")
			if target_species=="In" or target_species=="Ga":
				if pris_atoms[index]=="P":
					for j,third_layer_index in enumerate(pris_third_layer[count]):
						out.write(str(pris_atoms[third_layer_index]))
						out.write("  ")
						out.write(str(pris_coords[third_layer_index][0]+.392))
						out.write("  ")
						out.write(str(pris_coords[third_layer_index][1]))
						out.write("  ")
						out.write(str(pris_coords[third_layer_index][2]-.392))
						out.write(" \n")
					count=count+1
		else:
			out.write(str(pris_atoms[index]))
			out.write("  ")
			out.write(str(pris_coords[index][0]))
			out.write("  ")
			out.write(str(pris_coords[index][1]))
			out.write("  ")
			out.write(str(pris_coords[index][2]))
			out.write(" \n")
			if target_species=="In" or target_species=="Ga":
				if pris_atoms[index]=="P":
					for j,third_layer_index in enumerate(pris_third_layer[count]):
						out.write(str(pris_atoms[third_layer_index]))
						out.write("  ")
						out.write(str(pris_coords[third_layer_index][0]))
						out.write("  ")
						out.write(str(pris_coords[third_layer_index][1]))
						out.write("  ")
						out.write(str(pris_coords[third_layer_index][2]))
						out.write(" \n")
					count=count+1

	out.write(str(n_atoms)+"\n")
	out.write("Defect \n")
	for i,index in enumerate(def_targets):
		out.write(target_species)
		out.write("  ")
		out.write(str(def_coords[index][0]))
		out.write("  ")
		out.write(str(def_coords[index][1]))
		out.write("  ")
		out.write(str(def_coords[index][2]))
		out.write(" \n")
	count=0
	for i,index in enumerate(def_bonded):
		out.write(str(def_atoms[index]))
		out.write("  ")
		out.write(str(def_coords[index][0]))
		out.write("  ")
		out.write(str(def_coords[index][1]))
		out.write("  ")
		out.write(str(def_coords[index][2]))
		out.write(" \n")
		if target_species=="In" or target_species=="Ga":
			if def_atoms[index]=="P":
				for i,index2 in enumerate(third_layer[count]):
					out.write(str(def_atoms[index2]))
					out.write("  ")
					out.write(str(def_coords[index2][0]))
					out.write("  ")
					out.write(str(def_coords[index2][1]))
					out.write("  ")
					out.write(str(def_coords[index2][2]))
					out.write(" \n")
				count=count+1