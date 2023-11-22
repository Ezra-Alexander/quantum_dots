import numpy as np
import sys
import math
import random
from qchem_helper import read_xyz, write_xyz, get_geom_io, get_converged_geom
from geom_helper import get_nn, get_dists

#the sister script to the bash script of the same name
#inputs a file containing a geometry, an element to be replaced, an element to replace it with, and a number of replacements to make
#outputs a specified number of .xyz files, each of which has a random, unique version of the replacement

#relies on qchem_helper

#inputs
xyz_source=sys.argv[1]
orig=sys.argv[2]
replace=sys.argv[3]
nsub=int(sys.argv[4])

lig_cutoff=2.6 #a manual parameter for determining how close a ligand needs to be to a cation to be bonded to it

#get the original .xyz
#works with any kind of input (ideally) - .xyz, plot_optzd.in, or opt.out
dot=xyz_source.index(".")
suffix=xyz_source[dot+1:]
if suffix=="xyz":
	coords,atoms=read_xyz(xyz_source)
elif suffix=="in":
	atoms,coords=get_geom_io(xyz_source)
elif suffix=="out":
	coords,atoms=get_converged_geom(xyz_source)


#groundwork for substitution
orig_ind = atoms==orig
in_ind = atoms=="In" 
ga_ind = atoms=="Ga"
cat_ind=np.logical_or(in_ind,ga_ind)
an_ind = atoms=="P"
lig_ind = np.logical_not(np.logical_or(cat_ind,an_ind))
norig=np.count_nonzero(orig_ind)
sub_ind=[]
while len(sub_ind)<nsub:
	new_ind=random.randint(1,norig)
	if new_ind not in sub_ind:
		sub_ind.append(new_ind)

print("Replacing",orig,"atoms: ", sub_ind)

#get distances, nearest neighbors
dists,inp_dists_all,in_lig_dists_all, in_p_lig_dists_all, pin_dists_all=get_dists(coords,cat_ind,an_ind,lig_ind) #despite naming, works for systems other than InP
#lexie's code doesn't find nn for ligands ...
#need to do ligand nn myself



#the for loop that does the substitutions and writes the .xyzs
for i,ind in enumerate(sub_ind):
	orig_count=0
	new_atoms=[]
	new_xyz=[]
	for j,atom in enumerate(atoms):
		if orig_ind[j]:
			orig_count=orig_count+1

			if orig_count==ind:
				#do replacement

				#replace atomic ligand with hydroxide. Cation - O bond length is not changed, O-H bond length of 0.97 A is used, and a rough Cation - O - H bond angle of 108 degrees is used
				if (replace == "OH" or replace == "oh") and (orig=="F" or orig=="Cl" or orig == "Br" or orig=="H"):

					num_nn=np.sum(dists[j] < lig_cutoff)

					if num_nn==2: #terminal hydroxide (counts itself)
						print(str(i+1)+"th replacement is a Terminal Hydroxide")
						new_atoms.append("O")
						new_xyz.append(coords[j])

						#we need to find the cation closest to this one
						orig_dists=np.delete(dists[j],np.where(dists[j]==np.amin(dists[j]))[0])
						neighbor=coords[np.where(dists[j]==np.amin(orig_dists))[0]][0]
						
						vector=coords[j]-neighbor
						unit=vector/np.linalg.norm(vector)
						normal_unit = np.cross(unit,[1, 0, 0]) / np.linalg.norm(np.cross(unit,[1, 0, 0])) #sort of arbitrary, hopefully it works

						angle=math.radians(108-90)
						v_parallel=math.sin(angle)*0.97*unit
						v_perp=math.cos(angle)*0.97*normal_unit

						h_coords=coords[j]+v_parallel+v_perp
						new_atoms.append("H")
						new_xyz.append(h_coords)


					elif num_nn==3: #bridging hydroxide
						print(str(i+1)+"th replacement is a Bridging Hydroxide")
						new_atoms.append("O")
						new_xyz.append(coords[j])

						#now we need the two closest
						argsort=np.argsort(dists[j])
						neighbor1=coords[argsort[1]]
						neighbor2=coords[argsort[2]]

						angle=math.radians(108-90)
						v1=coords[j]-neighbor1
						v2=coords[j]-neighbor2
						v_parallel=v1+v2
						u_parallel=v_parallel/np.linalg.norm(v_parallel)
						#v_perp=np.cross(v1,v2)
						#u_perp=v_perp/np.linalg.norm(v_perp)

						h_coords=coords[j]+(0.97*u_parallel)
						new_atoms.append("H")
						new_xyz.append(h_coords)


					else:
						print("Something went wrong - change cutoff value?")

				else:
					print("Performing generic atomic substitution. If this wasn't what you wanted, either your replacement isn't supported yet or it wasn't formatted correctly.")
					new_atoms.append(replace)
					new_xyz.append(coords[j])

			else:
				new_atoms.append(atom)
				new_xyz.append(coords[j])

		else:
			new_atoms.append(atom)
			new_xyz.append(coords[j])

	write_xyz(str(i+1)+".xyz",new_atoms,new_xyz)