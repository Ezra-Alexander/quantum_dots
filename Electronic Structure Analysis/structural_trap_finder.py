import numpy as np
import sys
import matplotlib.pyplot as plt
import numpy.linalg as npl
import copy
from geom_helper import *
from pdos_helper import dos_grid_general,get_alpha,get_ao_ind, get_ind_ao_underc
from openpyxl import Workbook, load_workbook
import os

# the goal of this script is to automatically identify the structural traps in a given QD
# this is hard to do without orbital localization, so this is not meant to be perfect
# but the basic idea is to:
#	1. scan through the DOS for trap states that arise from "4c" atoms
#	2. catalog structural defects, including 1 or 2 centers
#	3. cutout .xyzs for each defect
# I'm just going to read in band edges, not dealing with that can of worms
# Meant to be as generalizable as possible
# This excludes the most delocalized "traps" from consideration
# Data is written to a local excel file. 
# A second script can export this excel data to a bigger file.
#		For ML purposes, you really need data on states that are not structural traps. But that can wait
# Assumes all ligands are -1 charge

xyz_file=sys.argv[1] #your .xyz file
bas_file = sys.argv[2]      # txt file with total number of orbitals and occupied orbitals
coeff_file = sys.argv[3]    # txt version of qchem 53.0 OR numpy version
ipr_file = sys.argv[4]	# csv file with the ipr (now formally Paticipation Ratio) for each MO 
occ_min=int(sys.argv[5]) #the relative index of the 1st occupied bulk state. If HOMO-17, for example, put 17
vir_max=int(sys.argv[6]) #the relative index of the 1st virtual bulk state. If LUMO+17, for example, put 17

alpha_thresh=0.9 #an arbitrary parameter (between 0 and 1) for determining what percent of the max individual contribution to an MO an atom must have before it can be considered trapping. A value of 1 would exclude all atoms that are not the the maximum contributor
double_thresh=0.666666 #an arbitrary parameter (between 0 and 1) for determining if a structural trap is one or two center. The lower the value, the more traps will be considered two-centered
ipr_thresh=0.4 #an arbitrary parameter (between 0 and 1) for determining if a state is "localized enough" to be a trap. The lower the value, the fewer states will be considered traps 

orb_per_atom_def2svp={'In': 26, 'Ga': 32, 'P': 18, 'Cl': 18, 'Br': 32,'F':14,'H':5,'O':14, 'C':14,'S':18, 'Li':9, 'Al':18, "Zn":31,"S":18,"Se":32, "Si":18}
orb_per_atom=orb_per_atom_def2svp # choose which dictionary to use
covalent_radii={"In":1.42,"Ga":1.22,"P":1.07,"Cl":1.02,"F":0.57,"Zn":1.22,"S":1.05,"Se":1.20,"O":0.66,"Al":1.21,"H":0.31}

dir_name="structural_cutouts" #the name of the directory the .xyz cutouts will be written to

# parse orbital info
nbas,nocc=np.loadtxt(bas_file,dtype=int,unpack=True)
homo=nocc-1

# read xyz file
coords,atoms=read_input_xyz(xyz_file)

# read 53.npz
coeff_file_expand=np.load(coeff_file)
mo_mat=coeff_file_expand['arr_0'] # normalized
mo_e = coeff_file_expand['arr_1']
mo_e = mo_e * 27.2114 # MO energy, in eV

#read low_orb_ipr.csv
ipr=np.loadtxt(ipr_file)

#step 1 is to identify the under-coordinated atoms in our QD
#the end result needs to be those boolean numpy arrays to feed into get_ind_ao_underc to feed into get_alpha

ind_In = (atoms=='In')
ind_P = np.logical_or((atoms=='P'), np.logical_or((atoms=='S'),(atoms=='Se')))
ind_Ga = (atoms=='Ga')
ind_Al=np.logical_or((atoms=="Al"),(atoms=="Zn"))
ind_cat = np.logical_or(ind_In, np.logical_or(ind_Ga,ind_Al))
ind_lig = np.logical_or(atoms=="F",np.logical_or((atoms=="Cl"),(atoms=="H")))
ind_other = np.logical_not(np.logical_or(ind_cat,np.logical_or(ind_P,ind_lig)))
if np.count_nonzero(ind_other)>0:
	raise Exception("There's an element here that isn't supported yet!")

dists=dist_all_points(coords)
connectivity=smart_connectivity_finder(dists,atoms)
angles=get_all_angles(coords)

ind_cat_uc=[]
ind_p_uc=[]
for i,atom in enumerate(atoms):
	if ind_cat[i] and len(connectivity[i])<4:
		ind_cat_uc.append(True)
		ind_p_uc.append(False)
	elif ind_P[i] and len(connectivity[i])<4:
		ind_p_uc.append(True)
		ind_cat_uc.append(False)
	else: 
		ind_cat_uc.append(False)
		ind_p_uc.append(False)
ind_cat_uc=np.array(ind_cat_uc)
ind_p_uc=np.array(ind_p_uc)

ind_cat_ao,ind_p_ao,ind_lig_ao = get_ao_ind([ind_cat,ind_P,ind_lig],atoms,nbas,orb_per_atom)
cat_underc_ind_ao,p_underc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,ind_cat_uc,ind_p_uc)

p_nouc_ao = np.logical_xor(p_underc_ind_ao,ind_p_ao)
cat_nouc_ao = np.logical_xor(cat_underc_ind_ao,ind_cat_ao)

alpha_cat_uc,alpha_cat_fc,alpha_p_uc,alpha_p_fc,alpha_lig=get_alpha(mo_mat,[cat_underc_ind_ao,cat_nouc_ao,p_underc_ind_ao,p_nouc_ao,ind_lig_ao])

test = np.all(np.isclose(alpha_cat_uc+alpha_cat_fc+alpha_p_uc+alpha_p_fc+alpha_lig,1))
if test == False: raise ValueError('Alpha doesnt add to 1!')
print('Alphas add to 1?:',test)

# I also need to hold on to information about what particular atoms contribute to each MO
# this is such a stupid way of doing this
every_atom=[]
for i,atom in enumerate(atoms):
	ind_i=np.full(len(atoms),False)
	ind_i[i]=True
	every_atom.append(ind_i)
every_atom=np.array(every_atom)

every_atom_ao=get_ao_ind(every_atom,atoms,nbas,orb_per_atom)
alpha_total=get_alpha(mo_mat,every_atom_ao)
atoms_by_orbs=np.array(alpha_total)
orbs_by_atoms=atoms_by_orbs.T

#now i want to loop over the orbitals, occupied then virtual, and use the alphas to identify 4c traps
#we want specifically traps where the main individual contribution comes from a 4c species
#but because of state mixing we may need a bit of a grace window
#I also want to exclude the most delocalized states
n_p=np.count_nonzero(ind_P)
n_cat=np.count_nonzero(ind_cat)
alpha_p=alpha_p_fc+alpha_p_uc
alpha_cat=alpha_cat_fc+alpha_cat_uc

occ_4c_traps=[]
occ_4c_centers=[]
orb=0
while orb<occ_min:
	mo_i=homo-orb

	ipr_i=ipr[mo_i]
	alpha_p_i=alpha_p[mo_i]
	nfrac=n_p/len(atoms)
	ipr_cutoff=(ipr_thresh*nfrac)/alpha_p_i
	if ipr_i<ipr_cutoff:

		mo_alphas=orbs_by_atoms[mo_i]
		max_alpha=np.max(mo_alphas)
		thresh=max_alpha*alpha_thresh
		centers=[]
		for i,atom in enumerate(atoms):
			if mo_alphas[i]>=thresh:
				if len(connectivity[i])>3:
					if mo_i not in occ_4c_traps:
						occ_4c_traps.append(mo_i)

					#need to prevent adding centers that are adjacent as "separate" centers
					test=True
					for j,neighbor in enumerate(connectivity[i]):
						if [neighbor] in centers:
							test=False
						for k,neighbor2 in enumerate(connectivity[neighbor]):
							if [neighbor2] in centers:
								test=False

					if test:
						centers.append([i])


		if len(centers)>0:
			occ_4c_centers.append(centers)
	orb=orb+1

print("Occupied structural traps:",[x+1 for x in occ_4c_traps])

vir_4c_traps=[]
vir_4c_centers=[]
orb=0
while orb<vir_max:
	mo_i=nocc+orb

	ipr_i=ipr[mo_i]
	alpha_cat_i=alpha_cat[mo_i]
	nfrac=n_cat/len(atoms)
	ipr_cutoff=(ipr_thresh*nfrac)/alpha_cat_i
	if ipr_i<ipr_cutoff:

		mo_alphas=orbs_by_atoms[mo_i]
		max_alpha=np.max(mo_alphas)
		thresh=max_alpha*alpha_thresh
		centers=[]
		for i,atom in enumerate(atoms):
			if mo_alphas[i]>=thresh:
				if len(connectivity[i])>3:
					if mo_i not in vir_4c_traps:

						vir_4c_traps.append(mo_i)

					#need to prevent adding centers that are adjacent as "separate" centers
					test=True
					for j,neighbor in enumerate(connectivity[i]):
						if [neighbor] in centers:
							test=False
						for k,neighbor2 in enumerate(connectivity[neighbor]):
							if [neighbor2] in centers:
								test=False
					if test:
						centers.append([i])

		if len(centers)>0:
			vir_4c_centers.append(centers)			
	orb=orb+1

print("Virtual structural traps:",[x+1 for x in vir_4c_traps])

structural_traps=occ_4c_traps+vir_4c_traps
structural_trap_centers=occ_4c_centers+vir_4c_centers

#Now, for each structural trap, we want to identify the structural motif that contributes to the trap
#We already have the centers
#I should write the .xyzs here

if not os.path.exists(dir_name):
  os.mkdir(dir_name)

cutouts=[]
for i,trap in enumerate(structural_traps):
	mo_alphas=orbs_by_atoms[trap]
	second_center=0
	for t, center in enumerate(structural_trap_centers[i]):
		center_i=center[0]
		#get coordination sphere
		next_nearest=[]
		for j,neighbor in enumerate(connectivity[center_i]):
			for k,neighbors in enumerate(connectivity[neighbor]):
				if (neighbors not in next_nearest) and (neighbors!=center_i):
					next_nearest.append(neighbors)
		full_sphere=connectivity[center_i]+next_nearest

		#check if any members of coordination sphere contribute similarly to center
		thresh=mo_alphas[center_i]*double_thresh
		for j,neighbor in enumerate(full_sphere):
			if mo_alphas[neighbor]>=thresh and len(connectivity[neighbor])>=4:
				print("Trap",trap+1,"is a two center trap with contributions from both",atoms[center_i],to_atom_specific_index(atoms,center_i),"and",atoms[neighbor],to_atom_specific_index(atoms,neighbor))
				second_center=second_center+1
				if second_center==1:
					center_2_i=neighbor
					structural_trap_centers[i][t].append(neighbor)
				else:
					raise Exception("You have a (currently unsupported) 3-center trap in state",trap+1, "on atoms:",atoms[center_i],to_atom_specific_index(atoms,center_i),atoms[center_2_i],to_atom_specific_index(atoms,center_2_i),atoms[neighbor],to_atom_specific_index(atoms,neighbor))


		#write centers and their nearest neighbors to a .xyz
		new_coords=[]
		new_atoms=[]
		new_coords.append(coords[center_i])
		new_atoms.append(atoms[center_i])
		if second_center==1:
			new_coords.append(coords[center_2_i])
			new_atoms.append(atoms[center_2_i])
			for j,neighbor in enumerate(connectivity[center_i]):
				if neighbor!=center_2_i:
					new_coords.append(coords[neighbor])
					new_atoms.append(atoms[neighbor])
			for j,neighbor in enumerate(connectivity[center_2_i]):
				if neighbor not in connectivity[center_i] and neighbor!=center_i:
					new_coords.append(coords[neighbor])
					new_atoms.append(atoms[neighbor])
		elif second_center==0:
			for j,neighbor in enumerate(connectivity[center_i]):
				new_coords.append(coords[neighbor])
				new_atoms.append(atoms[neighbor])

		write_xyz(dir_name+"/"+str(trap+1)+"_cutout_"+str(t+1)+".xyz",new_atoms,new_coords)
		cutouts.append([new_atoms,new_coords])

#Finally, let's write all the relevant info to an excel file

wb=Workbook()
ws1 = wb.active #Just the one worksheet

#write column labels
#there's no general dot information here - I can add that when I combine this with the mothership
#this format supports atoms up to 6-coordinate
#Booleans should be set as 0 or 1

ws1["A1"]="MO Index" #presumably to be excluded from training
ws1["B1"]="Occupied?"
ws1["C1"]="Total energy (eV)"
ws1["D1"]="Depth (eV)"
ws1["E1"]="Participation Ratio"
ws1["F1"]="1 or 2 centers?"
ws1["G1"]="Inter-Center Distance (A)"
ws1["H1"]="Total # of Atoms"
ws1["I1"]="Center 1 - Element"
ws1["J1"]="Center 1 - Element-Specific Index"
ws1["K1"]="Center 1 - Overall Index"
ws1["L1"]="Center 1 - Coordination #"
ws1["M1"]="Center 1 - Alpha"
ws1["N1"]="Center 1 - X-Coordinate (A)"
ws1["O1"]="Center 1 - Y-Coordinate (A)"
ws1["P1"]="Center 1 - Z-Coordinate (A)"
ws1["Q1"]="Center 1 - Bonded Element 1"
ws1["R1"]="Center 1 - Bonded Element 2"
ws1["S1"]="Center 1 - Bonded Element 3"
ws1["T1"]="Center 1 - Bonded Element 4"
ws1["U1"]="Center 1 - Bonded Element 5"
ws1["V1"]="Center 1 - Bonded Element 6"
ws1["W1"]="Center 1 - Bonded Element 1 Under-Coordinated?" 
ws1["X1"]="Center 1 - Bonded Element 2 Under-Coordinated?"
ws1["Y1"]="Center 1 - Bonded Element 3 Under-Coordinated?"
ws1["Z1"]="Center 1 - Bonded Element 4 Under-Coordinated?"
ws1["AA1"]="Center 1 - Bonded Element 5 Under-Coordinated?"
ws1["AB1"]="Center 1 - Bonded Element 6 Under-Coordinated?"
ws1["AC1"]="Center 1 - Bond Length 1 (A)"
ws1["AD1"]="Center 1 - Bond Length 2 (A)"
ws1["AE1"]="Center 1 - Bond Length 3 (A)"
ws1["AF1"]="Center 1 - Bond Length 4 (A)"
ws1["AG1"]="Center 1 - Bond Length 5 (A)"
ws1["AH1"]="Center 1 - Bond Length 6 (A)"
ws1["AI1"]="Center 1 - Bond Length 1 Ratio to Average" #"average" is just sum of covalent radii
ws1["AJ1"]="Center 1 - Bond Length 2 Ratio to Average"
ws1["AK1"]="Center 1 - Bond Length 3 Ratio to Average"
ws1["AL1"]="Center 1 - Bond Length 4 Ratio to Average"
ws1["AM1"]="Center 1 - Bond Length 5 Ratio to Average"
ws1["AN1"]="Center 1 - Bond Length 6 Ratio to Average"
ws1["AO1"]="Center 1 - Bond Angle 1 (Elements 1,2)"
ws1["AP1"]="Center 1 - Bond Angle 2 (Elements 1,3)"
ws1["AQ1"]="Center 1 - Bond Angle 3 (Elements 1,4)"
ws1["AR1"]="Center 1 - Bond Angle 4 (Elements 1,5)"
ws1["AS1"]="Center 1 - Bond Angle 5 (Elements 1,6)"
ws1["AT1"]="Center 1 - Bond Angle 6 (Elements 2,3)"
ws1["AU1"]="Center 1 - Bond Angle 7 (Elements 2,4)"
ws1["AV1"]="Center 1 - Bond Angle 8 (Elements 2,5)"
ws1["AW1"]="Center 1 - Bond Angle 9 (Elements 2,6)"
ws1["AX1"]="Center 1 - Bond Angle 10 (Elements 3,4)"
ws1["AY1"]="Center 1 - Bond Angle 11 (Elements 3,5)"
ws1["AZ1"]="Center 1 - Bond Angle 12 (Elements 3,6)"
ws1["BA1"]="Center 1 - Bond Angle 13 (Elements 4,5)"
ws1["BB1"]="Center 1 - Bond Angle 14 (Elements 4,6)"
ws1["BC1"]="Center 1 - Bond Angle 15 (Elements 5,6)"
ws1["BD1"]="Center 2 - Element"
ws1["BE1"]="Center 2 - Element-Specific Index"
ws1["BF1"]="Center 2 - Overall Index"
ws1["BG1"]="Center 2 - Coordination #"
ws1["BH1"]="Center 2 - Alpha"
ws1["BI1"]="Center 2 - X-Coordinate (A)"
ws1["BJ1"]="Center 2 - Y-Coordinate (A)"
ws1["BK1"]="Center 2 - Z-Coordinate (A)"
ws1["BL1"]="Center 2 - Bonded Element 1"
ws1["BM1"]="Center 2 - Bonded Element 2"
ws1["BN1"]="Center 2 - Bonded Element 3"
ws1["BO1"]="Center 2 - Bonded Element 4"
ws1["BP1"]="Center 2 - Bonded Element 5"
ws1["BQ1"]="Center 2 - Bonded Element 6"
ws1["BR1"]="Center 2 - Bonded Element 1 Under-Coordinated?" 
ws1["BS1"]="Center 2 - Bonded Element 2 Under-Coordinated?"
ws1["BT1"]="Center 2 - Bonded Element 3 Under-Coordinated?"
ws1["BU1"]="Center 2 - Bonded Element 4 Under-Coordinated?"
ws1["BV1"]="Center 2 - Bonded Element 5 Under-Coordinated?"
ws1["BW1"]="Center 2 - Bonded Element 6 Under-Coordinated?"
ws1["BX1"]="Center 2 - Bond Length 1 (A)"
ws1["BY1"]="Center 2 - Bond Length 2 (A)"
ws1["BZ1"]="Center 2 - Bond Length 3 (A)"
ws1["CA1"]="Center 2 - Bond Length 4 (A)"
ws1["CB1"]="Center 2 - Bond Length 5 (A)"
ws1["CC1"]="Center 2 - Bond Length 6 (A)"
ws1["CD1"]="Center 2 - Bond Length 1 Ratio to Average"
ws1["CE1"]="Center 2 - Bond Length 2 Ratio to Average"
ws1["CF1"]="Center 2 - Bond Length 3 Ratio to Average"
ws1["CG1"]="Center 2 - Bond Length 4 Ratio to Average"
ws1["CH1"]="Center 2 - Bond Length 5 Ratio to Average"
ws1["CI1"]="Center 2 - Bond Length 6 Ratio to Average"
ws1["CJ1"]="Center 2 - Bond Angle 1 (Elements 1,2)"
ws1["CK1"]="Center 2 - Bond Angle 2 (Elements 1,3)"
ws1["CL1"]="Center 2 - Bond Angle 3 (Elements 1,4)"
ws1["CM1"]="Center 2 - Bond Angle 4 (Elements 1,5)"
ws1["CN1"]="Center 2 - Bond Angle 5 (Elements 1,6)"
ws1["CO1"]="Center 2 - Bond Angle 6 (Elements 2,3)"
ws1["CP1"]="Center 2 - Bond Angle 7 (Elements 2,4)"
ws1["CQ1"]="Center 2 - Bond Angle 8 (Elements 2,5)"
ws1["CR1"]="Center 2 - Bond Angle 9 (Elements 2,6)"
ws1["CS1"]="Center 2 - Bond Angle 10 (Elements 3,4)"
ws1["CT1"]="Center 2 - Bond Angle 11 (Elements 3,5)"
ws1["CU1"]="Center 2 - Bond Angle 12 (Elements 3,6)"
ws1["CV1"]="Center 2 - Bond Angle 13 (Elements 4,5)"
ws1["CW1"]="Center 2 - Bond Angle 14 (Elements 4,6)"
ws1["CX1"]="Center 2 - Bond Angle 15 (Elements 5,6)"
ws1["CY1"]="Employed Max Alpha Contribution Threshold"
ws1["CZ1"]="Employed Alpha 2-Center Threshold"
ws1["DA1"]="Employed State-Specific IPR Threshold"
ws1["DB1"]="Mixed with an Under-Coordinated Trap?" #is uc density / fc density greater than the atomic ratio?
ws1["DC1"]="Total Under-Coordinated Cation Alpha"
ws1["DD1"]="Total Fully-Coordinated Cation Alpha"
ws1["DE1"]="Total Under-Coordinated Anion Alpha"
ws1["DF1"]="Total Fully-Coordinated Anion Alpha"
ws1["DG1"]="Total Ligand Alpha"
#dot-specific info. I'll append the name to the front
ws1["DH1"]="Total # Cations"
ws1["DI1"]="Total # Anions"
ws1["DJ1"]="Total # Ligand"
ws1["DK1"]="Fraction of Cations that are In"
ws1["DL1"]="Total Charge"
ws1["DM1"]="HOMO-LUMO Gap (eV)"
ws1["DN1"]="VBM relative index"
ws1["DO1"]="CBM relative index"
ws1["DP1"]="Bulk-Bulk Band Gap (eV)"
ws1["DQ1"]="Quantum Confined State?"




#now let's actually fill in the rows
row="1"
cutout_count=0
for i,trap in enumerate(structural_traps):
	for j,defect in enumerate(structural_trap_centers[i]):
		row=str(int(row)+1)

		ws1["A"+row]=trap+1

		if trap<=homo:
			ws1["B"+row]=1
			ws1["D"+row]=mo_e[trap]-mo_e[homo-occ_min]
		else:
			ws1["B"+row]=0
			ws1["D"+row]=mo_e[homo+1+vir_max]-mo_e[trap]

		ws1["C"+row]=mo_e[trap]
		ws1["E"+row]=ipr[trap]
		ws1["F"+row]=len(defect)

		if len(defect)>1:
			ws1["G"+row]= dists[defect[0]][defect[1]]

		defect_atoms=cutouts[cutout_count][0]
		defect_coords=cutouts[cutout_count][1]
		ws1["H"+row]= len(defect_atoms)
		cutout_count=cutout_count+1

		ws1["I"+row] = atoms[defect[0]]
		ws1["J"+row] = to_atom_specific_index(atoms,defect[0])
		ws1["K"+row] = defect[0]+1
		ws1["L"+row] = len(connectivity[defect[0]])
		ws1["M"+row] = orbs_by_atoms[trap][defect[0]]

		ws1["N"+row] = coords[defect[0]][0]
		ws1["O"+row] = coords[defect[0]][1]
		ws1["P"+row] = coords[defect[0]][2]

		col="Q"
		for k,bond in enumerate(connectivity[defect[0]]):
			neighbor=atoms[bond]
			ws1[col+row] = neighbor
			col=chr(ord(col)+1)

		col="W"
		for k,bond in enumerate(connectivity[defect[0]]):
			if (not ind_lig[bond]) and len(connectivity[bond])<4:
				ws1[col+row] = 1
			else:
				ws1[col+row] = 0

			one_letter=True
			if col=="Z":
				col="AA"
				one_letter=False
			elif one_letter:
				col=chr(ord(col)+1)
			else:
				col=col[0]+chr(ord(col[1])+1)

		col="AC"
		for k,bond in enumerate(connectivity[defect[0]]):
			ws1[col+row]=dists[bond][defect[0]]

			col=col[0]+chr(ord(col[1])+1)

		col="AI"
		for k,bond in enumerate(connectivity[defect[0]]):
			ws1[col+row]=dists[bond][defect[0]]/( covalent_radii[ atoms[bond] ] + covalent_radii[ atoms[defect[0]] ] )

			col=col[0]+chr(ord(col[1])+1)

		col="AO"
		for k,bond in enumerate(connectivity[defect[0]]):
			for l,bond2 in enumerate(connectivity[defect[0]]):
				if l>k:
					ws1[col+row]=angles[defect[0]][bond][bond2]

					if l==len(connectivity[defect[0]])-1:
						if k==0:
							col="AT"
						elif k==1:
							col="AX"
						elif k==2:
							col="BA"
						elif k==3:
							col="BC"
					elif col[1]=="Z":
						col=chr(ord(col[0])+1)+"A"
					else:
						col=col[0]+chr(ord(col[1])+1)

		if len(defect)>1:
			ws1["BD"+row]=atoms[defect[1]]
			ws1["BE"+row] = to_atom_specific_index(atoms,defect[1])
			ws1["BF"+row] = defect[1]+1
			ws1["BG"+row] = len(connectivity[defect[1]])
			ws1["BH"+row] = orbs_by_atoms[trap][defect[1]]

			ws1["BI"+row] = coords[defect[1]][0]
			ws1["BJ"+row] = coords[defect[1]][1]
			ws1["BK"+row] = coords[defect[1]][2]

			col="BL"
			for k,bond in enumerate(connectivity[defect[1]]):
				neighbor=atoms[bond]
				ws1[col+row] = neighbor

				col=col[0]+chr(ord(col[1])+1)

			col="BR"
			for k,bond in enumerate(connectivity[defect[1]]):
				if (not ind_lig[bond]) and len(connectivity[bond])<4:
					ws1[col+row] = 1
				else:
					ws1[col+row] = 0

				col=col[0]+chr(ord(col[1])+1)
				

			col="BX"
			for k,bond in enumerate(connectivity[defect[1]]):
				ws1[col+row]=dists[bond][defect[1]]

				if col[1]=="Z":
						col=chr(ord(col[0])+1)+"A"
				else:
					col=col[0]+chr(ord(col[1])+1)

			col="CD"
			for k,bond in enumerate(connectivity[defect[1]]):
				ws1[col+row]=dists[bond][defect[1]]/( covalent_radii[ atoms[bond] ] + covalent_radii[ atoms[defect[1]] ] )

				col=col[0]+chr(ord(col[1])+1)

			col="CJ"
			for k,bond in enumerate(connectivity[defect[1]]):
				for l,bond2 in enumerate(connectivity[defect[1]]):
					if l>k:
						ws1[col+row]=angles[defect[1]][bond][bond2]

						if l==len(connectivity[defect[1]])-1:
							if k==0:
								col="CO"
							elif k==1:
								col="CS"
							elif k==2:
								col="CV"
							elif k==3:
								col="CX"
						elif col[1]=="Z":
							col=chr(ord(col[0])+1)+"A"
						else:
							col=col[0]+chr(ord(col[1])+1)

		ws1["CY"+row]=alpha_thresh
		ws1["CZ"+row]=double_thresh
		ws1["DA"+row]=ipr_thresh


		if trap<=homo:
			if (alpha_p_uc[trap]/alpha_p_fc[trap])>(np.count_nonzero(ind_p_uc)/(np.count_nonzero(ind_P)-np.count_nonzero(ind_p_uc))):
				ws1["DB"+row]=1
			else:
				ws1["DB"+row]=0
		else:
			if (alpha_cat_uc[trap]/alpha_cat_fc[trap])>(np.count_nonzero(ind_cat_uc)/(np.count_nonzero(ind_cat)-np.count_nonzero(ind_cat_uc))):
				ws1["DB"+row]=1
			else:
				ws1["DB"+row]=0

		ws1["DC"+row]=alpha_cat_uc[trap]
		ws1["DD"+row]=alpha_cat_fc[trap]
		ws1["DE"+row]=alpha_p_uc[trap]
		ws1["DF"+row]=alpha_p_fc[trap]
		ws1["DG"+row]=alpha_lig[trap]

		ws1["DH"+row]=n_cat
		ws1["DI"+row]=n_p
		ws1["DJ"+row]=np.count_nonzero(ind_lig)
		ws1["DK"+row]=np.count_nonzero(ind_In)/n_cat
		ws1["DL"+row]=(3*n_cat)-(3*n_p)-np.count_nonzero(ind_lig)
		ws1["DM"+row]=mo_e[homo+1]-mo_e[homo]
		ws1["DN"+row]=occ_min
		ws1["DO"+row]=vir_max
		ws1["DP"+row]=mo_e[homo+1+vir_max]-mo_e[homo-occ_min]

		ipr_cutoff=(ipr_thresh*(n_cat/len(atoms)))/alpha_cat[homo]
		if ipr[homo+1]>ipr_cutoff:
			ws1["DQ"+row]=1
		else:
			ws1["DQ"+row]=0


wb.save("structural_trap_data.xlsx")