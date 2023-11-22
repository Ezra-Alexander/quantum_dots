import numpy as np
import sys
import matplotlib.pyplot as plt
import numpy.linalg as npl
from qd_helper import *
import copy
from geom_helper import *
from pdos_helper import dos_grid_general,get_alpha,get_ao_ind, get_ind_ao_underc

#sort of like the dos componenter code, except it finds the contribution from a set of atoms (of the same type) over a range of MOs 
#note that it always prints orbitals lowest to highest

QD_file=sys.argv[1] #your .xyz file
bas_file = sys.argv[2]      # txt file with total number of orbitals and occupied orbitals
coeff_file = sys.argv[3]    # txt version of qchem 53.0 OR numpy version

#range_min=int(sys.argv[4]) #the lowest MO index to print
#range_max=int(sys.argv[5]) #the highest MO index to print
occ=sys.argv[4] #alternate inputs for band edges. o for occupied orbitals, u for unoccupied
bulk_index=int(sys.argv[5]) #the relative index of the 1st bulk state. If HOMO-17, for example, put 17

element = sys.argv[6] #no overall indices option we die like men #currently works only for In, Ga, and P targets
indices=[int(x) for x in sys.argv[7:]] #just a list of element specific indices (don't repeat the elements)

nbas,nocc=np.loadtxt("nbas.txt",dtype=int,unpack=True) #only runs if nbas.txt is there
if occ=="o" or occ=="O":
	bulk=nocc-bulk_index
elif occ=="u" or occ=="U":
	bulk=nocc+1+bulk_index
else:
	print("4th input invalid: should be 'o' for occupied or 'u' for unoccupied")
#if range_min < nocc:
#	bulk=range_min
#else:
#	bulk=range_max

coords,atoms = read_input_xyz(QD_file)

nbas,nocc=np.loadtxt(bas_file,dtype=int,unpack=True)
homo=nocc-1

orb_per_atom={'In': 26, 'Ga': 32, 'P': 18, 'Cl': 18, 'Br': 32,'F':14,'H':5,'O':14, 'C':14,'S':18, 'Li':9, 'Al':18,'Zn':31,'S':18,"Se":32, "Si":18}

coeff_file_expand=np.load(coeff_file)
mo_mat=coeff_file_expand['arr_0'] # normalized
mo_e = coeff_file_expand['arr_1']

ind_In = atoms=='In'
ind_P = np.logical_or(atoms=='P',np.logical_or(atoms=="S",atoms=="Se"))
ind_Ga = atoms=="Ga"
ind_Al = np.logical_or(atoms=="Zn", atoms=="Al")
ind_lig = np.logical_or(atoms == "F",atoms=="Cl")
ind_lig2 = np.logical_or(np.logical_or(atoms=="Si",atoms=="O"),atoms=="C")

ga = False
if np.any(ind_Ga):
	ga = True

al = False
if np.any(ind_Al):
	al = True

#math
ind_in_ao,ind_p_ao,ind_ga_ao,ind_al_ao,ind_lig_ao,ind_lig2_ao = get_ao_ind([ind_In,ind_P,ind_Ga,ind_Al,ind_lig,ind_lig2],atoms,nbas,orb_per_atom)
'''
get squared coefficients on core,shell
'''
alpha_list = get_alpha(mo_mat,[ind_in_ao,ind_p_ao,ind_ga_ao,ind_al_ao,ind_lig_ao,ind_lig2_ao])
alpha_in,alpha_p,alpha_ga,alpha_al,alpha_lig,alpha_lig2=alpha_list
alpha_inga=alpha_in+alpha_ga

# make sure alphas add to 1
test = np.all(np.isclose(alpha_in+alpha_p+alpha_ga+alpha_al+alpha_lig+alpha_lig2,1))
if test == False: raise ValueError('Alpha doesnt add to 1!')
print('Alphas add to 1?:',test)
'''
calculate projected DOS
'''
mo_e = mo_e * 27.2114 # MO energy, in eV
E_grid = np.arange(-50,50,0.01) # energy grid to evaluate the DOS over
sigma=0.1 # broadening parameter
print('HOMO energy:',round(mo_e[homo],3),"eV")
print('HOMO-LUMO gap:', round(mo_e[homo+1]-mo_e[homo],3), "eV")
dos_list = dos_grid_general(E_grid, sigma,mo_e, alpha_list)
in_dos,p_dos,ga_dos,al_dos,lig_dos,lig2_dos=dos_list
inga_dos=in_dos+ga_dos


#so far that was all the normal PDOS code
#now make the "undercoordinated" index arrays
in_index_uc=np.zeros(len(atoms),dtype=bool)
p_index_uc=np.zeros(len(atoms),dtype=bool)
if ga:
	ga_index_uc=np.zeros(len(atoms),dtype=bool)
if al:
	al_index_uc=np.zeros(len(atoms),dtype=bool)

count=0
for i,atom in enumerate(atoms):
	if atom==element:
		count=count+1
		if count in indices:
			if element=="In":
				in_index_uc[i]=True
			#elif element=="P" or element=="S" or element=="Se":
			elif element=="P" or element=="S" or element=="Se":
				p_index_uc[i]=True
			elif ind_Ga[i]:
				ga_index_uc[i]=True
			elif element=="Al" or element=="Zn":
				al_index_uc[i]=True
			else:
				print("Your specified element isn't supported yet")

#do the rest of the PDOS-ing
in_uc_ind_ao,p_uc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,in_index_uc,p_index_uc)
if ga:
	ga_uc_ind_ao, p_uc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,ga_index_uc,p_index_uc)
if al:
	al_uc_ind_ao, p_uc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,al_index_uc,p_index_uc)

p_underc = p_uc_ind_ao #np.logical_or(sshell_underc_ind_ao,sshell_underc_ind_amb_ao)
in_underc = in_uc_ind_ao #np.logical_or(cdshell_underc_ind_ao,cdshell_underc_ind_amb_ao)
if ga:
	ga_underc = ga_uc_ind_ao
if al:
	al_underc = al_uc_ind_ao

p_nouc = np.logical_xor(p_underc,ind_p_ao)
in_nouc = np.logical_xor(in_underc,ind_in_ao)
if ga:
	ga_nouc = np.logical_xor(ga_underc,ind_ga_ao)
if al:
	al_nouc = np.logical_xor(al_underc,ind_al_ao)

alpha_in_uc,alpha_p_uc,alpha_in_fc,alpha_p_fc = get_alpha(mo_mat,[in_underc,p_underc,in_nouc,p_nouc])
alpha_fc=alpha_in_fc+alpha_p_fc
if ga:
	alpha_in_uc,alpha_p_uc,alpha_ga_uc,alpha_in_fc,alpha_p_fc,alpha_ga_fc = get_alpha(mo_mat,[in_underc,p_underc,ga_underc,in_nouc,p_nouc,ga_nouc])
	alpha_fc=alpha_in_fc+alpha_p_fc+alpha_ga_fc
if al:
	alpha_in_uc,alpha_p_uc,alpha_ga_uc,alpha_al_uc,alpha_in_fc,alpha_p_fc,alpha_ga_fc,alpha_al_fc = get_alpha(mo_mat,[in_underc,p_underc,ga_underc,al_underc,in_nouc,p_nouc,ga_nouc,al_nouc])
	alpha_fc=alpha_in_fc+alpha_p_fc+alpha_ga_fc+alpha_al_fc


#now print what we came for
print("MO Index   Target DOS   Remaining",element,"DOS   Energy (eV)   Depth from Bulk (eV)")
#if element=="In":
#	for i in range(range_max+1-range_min):
#		print(range_min+i,"     ",round(float(alpha_in_uc[range_min-1+i]),3),"      ", round(float(alpha_in_fc[range_min-1+i]),3),"             ",round(float(mo_e[range_min-1+i]),3), "     ", round(abs(float(mo_e[range_min-1+i])-float(mo_e[bulk-1])),3))
#elif element=="P":
#	for i in range(range_max+1-range_min):
#		print(range_min+i,"     ",round(float(alpha_p_uc[range_min-1+i]),3), "      ",round(float(alpha_p_fc[range_min-1+i]),3),"            ",round(float(mo_e[range_min-1+i]),3),"     ", round(abs(float(mo_e[range_min-1+i])-float(mo_e[bulk-1])),3))
#elif element=="Ga": 
#	for i in range(range_max+1-range_min):
#		print(range_min+i,"     ",round(float(alpha_ga_uc[range_min-1+i]),3), "      ",round(float(alpha_ga_fc[range_min-1+i]),3),"            ",round(float(mo_e[range_min-1+i]),3),"     ", round(abs(float(mo_e[range_min-1+i])-float(mo_e[bulk-1])),3))
#elif element=="Al": 
#	for i in range(range_max+1-range_min):
#		print(range_min+i,"     ",round(float(alpha_al_uc[range_min-1+i]),3), "      ",round(float(alpha_al_fc[range_min-1+i]),3),"            ",round(float(mo_e[range_min-1+i]),3),"     ", round(abs(float(mo_e[range_min-1+i])-float(mo_e[bulk-1])),3))

if occ=="o" or occ=="O":
	if element=="In":
		for i in range(bulk_index+1):
			print(bulk+i,"     ",round(float(alpha_in_uc[bulk-1+i]),3),"      ", round(float(alpha_in_fc[bulk-1+i]),3),"             ",round(float(mo_e[bulk-1+i]),3), "     ", round(abs(float(mo_e[bulk-1+i])-float(mo_e[bulk-1])),3))
	elif element=="P" or element=="S" or element=="Se":
		for i in range(bulk_index+1):
			print(bulk+i,"     ",round(float(alpha_p_uc[bulk-1+i]),3), "      ",round(float(alpha_p_fc[bulk-1+i]),3),"            ",round(float(mo_e[bulk-1+i]),3),"     ", round(abs(float(mo_e[bulk-1+i])-float(mo_e[bulk-1])),3))
	elif element=="Ga": 
		for i in range(bulk_index+1):
			print(bulk+i,"     ",round(float(alpha_ga_uc[bulk-1+i]),3), "      ",round(float(alpha_ga_fc[bulk-1+i]),3),"            ",round(float(mo_e[bulk-1+i]),3),"     ", round(abs(float(mo_e[bulk-1+i])-float(mo_e[bulk-1])),3))
	elif element=="Al": 
		for i in range(bulk_index+1):
			print(bulk+i,"     ",round(float(alpha_al_uc[bulk-1+i]),3), "      ",round(float(alpha_al_fc[bulk-1+i]),3),"            ",round(float(mo_e[bulk-1+i]),3),"     ", round(abs(float(mo_e[bulk-1+i])-float(mo_e[bulk-1])),3))
elif occ=="u" or occ=="U":
	if element=="In":
		for i in range(bulk_index+1):
			print(nocc+1+i,"     ",round(float(alpha_in_uc[nocc+i]),3),"      ", round(float(alpha_in_fc[nocc+i]),3),"             ",round(float(mo_e[nocc+i]),3), "     ", round(abs(float(mo_e[nocc+i])-float(mo_e[bulk-1])),3))
	elif element=="P" or element=="S":
		for i in range(bulk_index+1):
			print(nocc+1+i,"     ",round(float(alpha_p_uc[nocc+i]),3), "      ",round(float(alpha_p_fc[nocc+i]),3),"            ",round(float(mo_e[nocc+i]),3),"     ", round(abs(float(mo_e[nocc+i])-float(mo_e[bulk-1])),3))
	elif element=="Ga": 
		for i in range(bulk_index+1):
			print(nocc+1+i,"     ",round(float(alpha_ga_uc[nocc+i]),3), "      ",round(float(alpha_ga_fc[nocc+i]),3),"            ",round(float(mo_e[nocc+i]),3),"     ", round(abs(float(mo_e[nocc+i])-float(mo_e[bulk-1])),3))
	elif element=="Al" or element=="Zn": 
		for i in range(bulk_index+1):
			print(nocc+1+i,"     ",round(float(alpha_al_uc[nocc+i]),3), "      ",round(float(alpha_al_fc[nocc+i]),3),"            ",round(float(mo_e[nocc+i]),3),"     ", round(abs(float(mo_e[nocc+i])-float(mo_e[bulk-1])),3))








