from geom_helper import get_dists, get_underc_index, get_nn
from qd_helper import read_input_xyz, write_xyz
import sys
import math
import numpy as np
from optparse import OptionParser

#this code should take a QD .xyz file and add a specified (atomic) ligand to every 3c-surface atom of a specified type until those cations are fully coordinated
#Fails if there are 2c or less atoms

#The orginal .xyz file that I am editing
xyz = sys.argv[1]

#The name of the .xyz file this will make
name = sys.argv[2]

#the atomic ligand to passivate with
lig=sys.argv[3]

#the surface atom to passivate
base=sys.argv[4]

#the bond distance between your ligand and cation
bond=float(sys.argv[5])

cutoff=3.0 #arbitary cutoff for two atoms to be considered bonded. May need to be changed for some systems

#read source .xyz
coords,atoms=read_input_xyz(xyz)

#find undercoordinated cations
ind_cat=[]
ind_an=[]
ind_lig=[]
for i,atom in enumerate(atoms):
	if atom =="In" or atom == "Ga" or atom=="Al":
		ind_cat.append(True)
		ind_an.append(False)
		ind_lig.append(False)
	elif atom=="P":
		ind_cat.append(False)
		ind_an.append(True)
		ind_lig.append(False)
	elif atom=="F" or atom=="Cl" or atom=="H":
		ind_lig.append(True)
		ind_cat.append(False)
		ind_an.append(False)
	else:
		print("There's an atom in your source .xyz that isn't supported!")

ind_cat = np.array(ind_cat)
ind_an = np.array(ind_an)
ind_lig = np.array(ind_lig)

all_dists,catan_dists,catlig_dists,catanlig_dists,ancat_dists = get_dists(coords,ind_cat,ind_an,ind_lig)

Natoms = len(ind_cat)
all_nn,cat_nn_anlig,an_nn_cat = get_nn(catanlig_dists,ancat_dists,ind_cat,ind_an,cutoff,Natoms,ind_lig)


new_atoms=np.copy(atoms)
new_coords=np.copy(coords)
for i,atom in enumerate(atoms):
	if atom==base:
		if all_nn[i]==3:
			#find the nearest neighbors of our 3c atom
			nearest=[]
			nearest_dists=[]
			max_index=1000
			for j,dist in enumerate(all_dists[i]):
				if i!=j:
					if len(nearest)<3:
						nearest.append(j)
						if len(nearest_dists)>0:
							if dist>max(nearest_dists):
								max_index=len(nearest)-1
						else:
							max_index=len(nearest)-1
						nearest_dists.append(dist)
					elif dist<max(nearest_dists):
						nearest.pop(max_index)
						nearest_dists.pop(max_index)
						if dist>max(nearest_dists):
							max_index=2
						else:
							if nearest_dists[0]>nearest_dists[1]:
								max_index=0
							else:
								max_index=1
						nearest.append(j)
						nearest_dists.append(dist)


			#vector		
			v1=coords[i]-coords[nearest[0]] 
			v2=coords[i]-coords[nearest[1]] 
			v3=coords[i]-coords[nearest[2]] 
			vsum=v1+v2+v3
			u=vsum/np.linalg.norm(vsum)

			adding_coords=coords[i]+(u*bond)

			new_atoms=np.append(new_atoms,lig)
			new_coords=np.append(new_coords,[adding_coords],axis=0)


			
						
write_xyz(name,new_atoms,new_coords)