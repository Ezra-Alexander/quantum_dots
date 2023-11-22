from geom_helper import get_dists, get_underc_index
from qd_helper import read_input_xyz, write_xyz
import sys
import math
import numpy as np
from optparse import OptionParser

#inputs
parser=OptionParser()

parser.add_option("-s","--source",dest="source", help="The .xyz file with the bulk crystal structure. Must be specified. Currently only In, Ga, and P atoms are supported, but others can be easily implemented when needed.", action="store",type="string")
#source = sys.argv[1] #source .xyz file for the bulk crystal

parser.add_option("-o","--output",dest="output", default="output.xyz", help="The name of the file you want the final QD saved to. Default is output.xyz",action="store",type="string")
#output = sys.argv[2] #name of output .xyz file for QD

parser.add_option("-c","--center",dest="center", default="com",help="What type of atom the QD should be centered on. Default is the center of mass, abbreviated com. Midpoint between Cat and An (abbreviated mid) is also supported. InP source should always use a non-COM centering",action="store",type="string")
#center = sys.argv[3] #how the dot should be centered. Either "com" or the element that it should be centered on. Note that the COM of the GaP bulk is a Gallium

parser.add_option("-p", "--shape", dest="shape", default="sphere", help="The desired shape of the output QD. Default is sphere, though tetrahedra are also supported",action="store",type="string")
#shape = sys.argv[4] #what shape the dot should have. Either "sphere" or "tetrahedral"

parser.add_option("-l","--ligand", dest="ligand", default="F", help="The atomic ligand to passivate the QD surface with. Default is F. Current cutoff for atoms to be considered bonded is 3.0 Angstroms",action="store",type="string")
#ligand = sys.argv[5] #the ligand to passivate the surface with

parser.add_option("-d","--diameter", dest="diameter",help="The target QD diameter. Must be specified. For tetrahedra, diameter is twice the distance between the COM and the center of each face",action="store",type="string")
#diameter = sys.argv[6] #rough diameter of target dot, in Angstroms

(options,args)=parser.parse_args()
source = options.source
output = options.output
center = options.center
shape = options.shape
ligand = options.ligand
diameter = options.diameter


#read source .xyz
coords,atoms=read_input_xyz(source)


#info the code has
masses = {"Ga":69.723, "In":114.82, "P":30.974}

#choose center of dot

#compute center of mass
cent =np.array([0.0])
mass = 0
for i,atom in enumerate(atoms):
	cent = cent + coords[i]*masses[atom]
	mass = mass + masses[atom]
	com = cent/mass

if center == "com":
	
	start = com 

elif center == "mid":

	min_dist_an = 1000
	min_dist_cat=1000
	closest_an=len(atoms)
	closest_cat=len(atoms)
	for i,atom in enumerate(atoms):
		dist=np.sqrt(np.sum((coords[i] - com)**2))
		if atom == "In" or atom == "Ga":
			if dist < min_dist_cat:
				min_dist_cat=dist
				closest_cat=i
		if atom == "P":
			if dist < min_dist_an:
				min_dist_an=dist
				closest_an=i

	if closest_cat ==len(atoms) or closest_an==len(atoms):
		print("Something went wrong with centering!")

	start = (1.5*coords[closest_an]+0.5*coords[closest_cat])/2

else: #i.e. if a specific atom type is chosen

	min_dist = 1000 #arbitrarily large numbers
	closest = len(atoms)

	for i,atom in enumerate(atoms):
		if atom == center:
			dist=np.sqrt(np.sum((coords[i] - com)**2))
			if dist < min_dist:
				closest = i
				min_dist=dist

	if closest ==len(atoms):
		print("That atom isn't in this structure, silly!")

	start = coords[closest]


#carve the dot of specified shape

#going to start with the default of spherical. Would like to ultimately implement Tetrahedral and maybe even a cube

carved_atoms=[]
carved_coords=[]
if shape == "sphere":

	for i,atom in enumerate(atoms):
		dist=np.sqrt(np.sum((coords[i] - start)**2))
		if dist<float(diameter)/2:
			carved_atoms.append(atom)
			carved_coords.append(coords[i])


elif shape[:3]=="tet":

	#Find index of atom closest to center
	min_dist=10000 #arbitrarily large numbers
	min_in = len(atoms) 
	for i,atom in enumerate(atoms):
		dist=np.sqrt(np.sum((coords[i] - start)**2))
		if dist < min_dist:
			min_dist=dist
			min_in=i

	#Find 4 atoms closest to that atom
	neighbor_dists=[]
	neighbor_indexes=[]
	for i,atom in enumerate(atoms):	
		dist=np.sqrt(np.sum((coords[i] - coords[min_in])**2))
		if len(neighbor_dists)<4 and i != min_in:
			neighbor_dists.append(dist)
			neighbor_indexes.append(i)
		elif i==min_in:
			print("Making a Tetrahedral Dot")
		elif dist < max(neighbor_dists) and len(neighbor_dists)==4 and i != min_in:
			neighbor_indexes.pop(neighbor_dists.index(max(neighbor_dists)))
			neighbor_dists.pop(neighbor_dists.index(max(neighbor_dists)))
			neighbor_dists.append(dist)
			neighbor_indexes.append(i)

	#make unit vectors
	u1=(coords[neighbor_indexes[0]]-coords[min_in])/(neighbor_dists[0])
	u2=(coords[neighbor_indexes[1]]-coords[min_in])/(neighbor_dists[1])
	u3=(coords[neighbor_indexes[2]]-coords[min_in])/(neighbor_dists[2])
	u4=(coords[neighbor_indexes[3]]-coords[min_in])/(neighbor_dists[3])

	#when the center isn't at (0,0,0), you also need the projections of the center to compare to
	r1=np.sum(u1*start) #so these are dot products - projections
	r2=np.sum(u2*start)
	r3=np.sum(u3*start)
	r4=np.sum(u4*start)


	#Trying the troy code approach ...
	for i,atom in enumerate(atoms):

		c1=np.sum(u1*coords[i]) #so these are dot products - projections
		c2=np.sum(u2*coords[i])
		c3=np.sum(u3*coords[i])
		c4=np.sum(u4*coords[i])

		if (c1 > r1-float(diameter)/3) and (c2 > r2-float(diameter)/3) and (c3 > r3-float(diameter)/3) and (c4 > r4-float(diameter)/3):
			carved_atoms.append(atom)
			carved_coords.append(coords[i])

else:
	print("That shape isn't supported :(")


#Now we prune any 1-coordinate anions (also 0 coordinate)

#indexes
ind_cat=[]
ind_an=[]
for i,atom in enumerate(carved_atoms):
	if atom =="In" or atom == "Ga":
		ind_cat.append(True)
		ind_an.append(False)
	elif atom=="P":
		ind_cat.append(False)
		ind_an.append(True)
	else:
		print("There's an atom in your source .xyz that isn't supported!")

ind_cat = np.array(ind_cat)
ind_an = np.array(ind_an)

#find 1-c atoms (just using an_1c)
blank=np.zeros_like(ind_cat,dtype=bool) #another artifact of adapting code written for QDs with ligands
cat_1c, an_1c = get_underc_index(carved_coords,ind_cat,ind_an,blank, blank, 3.0,2) #2nd to last variable is nearest neighbor distance, may need to be changed for certain systems

#remove the naughty ones
carved_atoms_2=[]
carved_coords_2=[]
n_an=0
for i,atom in enumerate(carved_atoms):
	if ind_an[i]:
		if not an_1c[n_an]:
			carved_atoms_2.append(atom)
			carved_coords_2.append(carved_coords[i])
		n_an=n_an+1
	elif ind_cat[i]:
		carved_atoms_2.append(atom)
		carved_coords_2.append(carved_coords[i])


#Now we repeat for 1-coordinate cations (also 0 coordinate)
#can't do the two together because you'll still end up with In/Ga-1c if you remove the the P-1c from an In, for example

#indexes
ind_cat=[]
ind_an=[]
for i,atom in enumerate(carved_atoms_2):
	if atom =="In" or atom == "Ga":
		ind_cat.append(True)
		ind_an.append(False)
	elif atom=="P":
		ind_cat.append(False)
		ind_an.append(True)
	else:
		print("There's an atom in your source .xyz that isn't supported!")

ind_cat = np.array(ind_cat)
ind_an = np.array(ind_an)

#find 1-c atoms (just using an_1c)
blank=np.zeros_like(ind_cat,dtype=bool) #another artifact of adapting code written for QDs with ligands
cat_1c, an_1c = get_underc_index(carved_coords_2,ind_cat,ind_an,blank, blank, 3.0,2) #2nd to last variable is nearest neighbor distance, may need to be changed for certain systems

#remove the naughty ones
pruned_atoms=[]
pruned_coords=[]
n_cat=0
for i,atom in enumerate(carved_atoms_2):
	if ind_an[i]:
		pruned_atoms.append(atom)
		pruned_coords.append(carved_coords_2[i])	
	elif ind_cat[i]:
		if not cat_1c[n_cat]:
			pruned_atoms.append(atom)
			pruned_coords.append(carved_coords_2[i])
		n_cat=n_cat+1


#now we turn 2-Coordinate anions into our ligand

#re-index for 2c
ind_cat=[]
ind_an=[]
for i,atom in enumerate(pruned_atoms):
	if atom =="In" or atom == "Ga":
		ind_cat.append(True)
		ind_an.append(False)
	elif atom=="P":
		ind_cat.append(False)
		ind_an.append(True)
	else:
		print("There's an atom in your source .xyz that isn't supported!")

ind_cat = np.array(ind_cat)
ind_an = np.array(ind_an)

#find 2c atoms
blank=np.zeros_like(ind_cat,dtype=bool) #another artifact of adapting code written for QDs with ligands
cat_2c, an_2c = get_underc_index(pruned_coords,ind_cat,ind_an,blank, blank, 3.0,3) #2nd to last variable is nearest neighbor distance, may need to be changed for certain systems

#convert 2c anions
conv_atoms=[]
conv_coords=[]
n_an=0
for i,atom in enumerate(pruned_atoms):
	if ind_an[i]:
		if an_2c[n_an]:
			conv_atoms.append(ligand)
			conv_coords.append(pruned_coords[i])
		else:
			conv_atoms.append(atom)
			conv_coords.append(pruned_coords[i])
		n_an=n_an+1
	else:
		conv_atoms.append(atom)
		conv_coords.append(pruned_coords[i])


#now we passivate any undercoordinated cation. We can still use the cat_2c index above for 2c.

#first, add one ligand to each 2c cation in a way that maintains the tetrahedral geometry as well as possible
#note that troy's code leaves In/Ga-3c, presumably by taking In/Ga-2c and only adding one F
#meaning that this section should be skipped to make dots in best accordance with the original
#yes, I know this is horibly inefficient, don't @ me

#n_cat=0
#to_append =[]
#for i,atom in enumerate(conv_atoms):
#	if ind_cat[i]:
#		if cat_2c[n_cat]:
#			#compute all distances, take two lowest
#			dists=[]
#			for j,atom in enumerate(conv_atoms):
#				if j!=i:
#					dist=np.sqrt(np.sum((conv_coords[j] - conv_coords[i])**2))
#					dists.append([dist,j])
#			dists =sorted(dists, key=lambda x: x[0])
#			p1=conv_coords[dists[0][1]]
#			p2=conv_coords[dists[1][1]]
#			cent=conv_coords[i]
#			min_dist=dists[0][0]
#			#compute tetrahedral 3rd position
#			v1=p1-cent
#			v2=p2-cent
#			n=np.cross(v1,v2)
#			n = n/np.linalg.norm(n)
#			p_mid = (p1+p2)/2
#			v_mid = cent-p_mid
#			v_mid=v_mid/np.linalg.norm(v_mid)
#			p_new = cent+(math.sin(math.radians(36.2644))*v_mid*min_dist)+(math.cos(math.radians(36.2644))*n*min_dist) #note that we use the cation-anion bond length here, as this structure will be optimized
#			if len(to_append)>0: #checking to make sure I don't add the same point twice
#				test = True
#				for j,current in enumerate(to_append):
#					if math.isclose(p_new[0],current[0],abs_tol=0.1) and math.isclose(p_new[1],current[1],abs_tol=0.1) and math.isclose(p_new[2],current[2],abs_tol=0.1):
#						test=False
#				if test:
#					to_append.append(p_new)
#			else:
#				to_append.append(p_new)
#		n_cat=n_cat+1
#
#for i,lig in enumerate(to_append):
#	conv_atoms.append(ligand)
#	conv_coords.append(lig)


#now we go through and passivate all three-coordinate cations
#if skipping the previous section, we also passivate In/Ga-2c to 3c here
#re-index for 3c
ind_cat=[]
ind_an=[]
ind_lig=[]
for i,atom in enumerate(conv_atoms):
	if atom =="In" or atom == "Ga":
		ind_cat.append(True)
		ind_an.append(False)
		ind_lig.append(False)
	elif atom=="P":
		ind_cat.append(False)
		ind_an.append(True)
		ind_lig.append(False)
	elif atom==ligand:
		ind_lig.append(True)
		ind_cat.append(False)
		ind_an.append(False)
	else:
		print("There's an atom in your source .xyz that isn't supported!")

ind_cat = np.array(ind_cat)
ind_an = np.array(ind_an)
ind_lig=np.array(ind_lig)

#find 3c atoms
cat_3c, an_3c = get_underc_index(conv_coords,ind_cat,ind_an,ind_lig, ind_lig, 3.0,4) #2nd to last variable is nearest neighbor distance, may need to be changed for certain systems

#do the thing
n_cat=0
append_2 =[]
for i,atom in enumerate(conv_atoms):
	if ind_cat[i]:
		if cat_3c[n_cat]:
			dists=[]
			for j,atom2 in enumerate(conv_atoms):
				if j!=i:
					dist=np.sqrt(np.sum((conv_coords[j] - conv_coords[i])**2))
					dists.append([dist,j])
			dists =sorted(dists, key=lambda x: x[0])
			p1=conv_coords[dists[0][1]]
			p2=conv_coords[dists[1][1]]
			p3=conv_coords[dists[2][1]]

			dummy = (p1+p2+p3)/3
			if dists[1][0]+0.5<dists[2][0]:
				dummy = (p1+p2)/2

			vector = conv_coords[i]-dummy
			mag = np.linalg.norm(vector)
			unit = vector/mag

			bond_length=0
			if atom=="Ga" and ligand=="F":
				bond_length=1.75
			elif atom=="In" and ligand=="F":
				bond_length=1.95
			elif atom=="In" and ligand=="Cl":
				bond_length=2.37
			else:
				print("This ligand-cation pair needs its bond length hard coded :(")

			p_new=conv_coords[i]+(unit*bond_length)

			if len(append_2)>0: #checking to make sure I don't add the same point twice
				test = True
				for j,current in enumerate(append_2):
					if np.array_equal(p_new,current):
						print("Something might be weird with this structure ... Double check")
						test=False
				if test:
					append_2.append(p_new)
			else:
				append_2.append(p_new)
		n_cat=n_cat+1

for i,lig in enumerate(append_2):
	conv_atoms.append(ligand)
	conv_coords.append(lig)


#provide output information
#re-index for the last time
ind_cat=[]
ind_an=[]
ind_lig=[]
cat_type=""
an_type=""
for i,atom in enumerate(conv_atoms):
	if atom =="In" or atom == "Ga":
		ind_cat.append(True)
		ind_an.append(False)
		ind_lig.append(False)
		cat_type=atom
	elif atom=="P":
		ind_cat.append(False)
		ind_an.append(True)
		ind_lig.append(False)
		an_type=atom
	elif atom==ligand:
		ind_lig.append(True)
		ind_cat.append(False)
		ind_an.append(False)
	else:
		print("There's an atom in your source .xyz that isn't supported!")

ind_cat = np.array(ind_cat)
ind_an = np.array(ind_an)
ind_lig=np.array(ind_lig)

print("This dot has", np.count_nonzero(ind_cat), cat_type, ",", np.count_nonzero(ind_an), an_type, ", and", np.count_nonzero(ind_lig), ligand)
print("This leaves it with a charge of ", 3*(np.count_nonzero(ind_cat)-np.count_nonzero(ind_an))-np.count_nonzero(ind_lig)) #assumes a III-V dot

write_xyz(output,conv_atoms,conv_coords)