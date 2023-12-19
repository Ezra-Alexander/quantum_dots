import numpy as np
import sys
import matplotlib.pyplot as plt
import numpy.linalg as npl
from qd_helper import *
import copy
from geom_helper import *
from pdos_helper import dos_grid_general,get_alpha,get_ao_ind, get_ind_ao_underc

#the goal of this code is to make a non-janky way of assessing the different atomic contributions to the DOS of each atom in a given MO
#nore that dos contributions from ligands are ignored
#the way this code was originally implemented was horribly slow. This version is much better

xyz_file=sys.argv[1] #your .xyz file
bas_file = sys.argv[2]      # txt file with total number of orbitals and occupied orbitals
coeff_file = sys.argv[3]    # txt version of qchem 53.0 OR numpy version
target = int(sys.argv[4])-1 #integer index of MO you want to decompose

orb_per_atom_def2svp={'In': 26, 'Ga': 32, 'P': 18, 'Cl': 18, 'Br': 32,'F':14,'H':5,'O':14, 'C':14,'S':18, 'Li':9, 'Al':18, "Zn":31,"S":18,"Se":32, "Si":18}
orb_per_atom=orb_per_atom_def2svp # choose which dictionary to use

# parse orbital info
nbas,nocc=np.loadtxt(bas_file,dtype=int,unpack=True)
homo=nocc-1

# read xyz file
coords,atoms=read_input_xyz(xyz_file)

# read 53.npz
coeff_file_expand=np.load(coeff_file)
mo_mat=coeff_file_expand['arr_0'] # normalized
mo_e = coeff_file_expand['arr_1']

ind_In = (atoms=='In')
ind_P = np.logical_or((atoms=='P'), np.logical_or((atoms=='S'),(atoms=='Se')))
ind_Ga = (atoms=='Ga')
ind_Al=np.logical_or((atoms=="Al"),(atoms=="Zn"))
ind_cat = np.logical_or(ind_In, np.logical_or(ind_Ga,ind_Al))
ind_lig = np.logical_or(atoms=="F",np.logical_or((atoms=="Cl"),(atoms=="H")))
ind_other = np.logical_not(np.logical_or(ind_cat,np.logical_or(ind_P,ind_lig)))
if np.count_nonzero(ind_other)>0:
	raise Exception("There's an element here that isn't supported yet!")

# I need information about what particular atoms contribute to each MO
# this is such a stupid way of doing this
every_atom=[]
for i,atom in enumerate(atoms):
	ind_i=np.full(len(atoms),False)
	ind_i[i]=True
	every_atom.append(ind_i)
every_atom=np.array(every_atom)

every_atom_ao=get_ao_ind(every_atom,atoms,nbas,orb_per_atom)
alpha_total=get_alpha(mo_mat,every_atom_ao)

for i,atom in enumerate(atoms):
	if alpha_total[i][target]>0.02:
		print("In MO",str(target),",",atom,to_atom_specific_index(atoms,i),"(total index",str(i+1)+")","has a DOS contribution of",round(float(alpha_total[i][target]),3))