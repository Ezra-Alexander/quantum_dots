from pyscf import gto, dft, lo
import numpy as np
from pyscf.lo.boys import Boys
from pyscf.lo.pipek import Pipek
from pyscf.tools import cubegen
import h5py
import sys
from pdos_helper import dos_grid_general,get_alpha,get_ao_ind, get_ind_ao_underc
from qd_helper import read_input_xyz, get_atom_specific_index

#a PySCF script that should read in 53.npz and perform an orbital localization of a specified type on a specified subset of orbitals
#It then writes cube files for the desired,localized MOs
#Maybe I can even use the MO coefficients to do the DOS analysis? A bit ambitious

#So, there should be no actual SCF run here. I could run this all on my laptop, ideally

#Of course, this is all dependent on PySCF being able to perform orbital localizations on systems with Indium

xyz=sys.argv[1] #the .xyz file
bas_file=sys.argv[2] #ya boy nbas.txt
coeff_file=sys.argv[3] #53.npz, from QChem and my normal manipulations
charge=int(sys.argv[4]) #yeah I guess this could theoretically be computed from the .xyz but this is easier
loc_type=sys.argv[5] #maybe I'm getting a bit ambitious, but it would be neat to be able to choose your localization method
loc_min=int(sys.argv[6]) #the minimum index of orbitals to localize (inclusive), relative to HOMO/LUMO.  For example, enter 11 for HOMO-11. No non-continuous sets this time.
loc_max=int(sys.argv[7]) #the max index of orbitals to localize (inclusive), relative to HOMO/LUMO. For example, enter 0 for HOMO. I could theoretically implement an energy window but I won't
#if loc_max is > than loc_min, then we're doing UMOS. If loc_min > loc_max, we're doing OMOs

nbas,nocc=np.loadtxt(bas_file,dtype=int,unpack=True)
homo=nocc-1

#imported PDOS stuff, to attempt DOS analysis of localized MOs
coords,atoms=read_input_xyz(xyz)
orb_per_atom={'In': 26, 'Ga': 32, 'P': 18, 'Cl': 18, 'Br': 32,'F':14,'H':5,'O':14, 'C':14,'S':18, 'Li':9}
ind_In = (atoms=='In')
ind_P = np.logical_or((atoms=='P'), (atoms=='S'))
ind_Ga = (atoms=='Ga')
ind_lig = np.logical_or((atoms=='F'), (atoms=='Cl'))
ind_lig2=np.logical_or((atoms=="C"),np.logical_or((atoms=="O"),(atoms=="H")))

ind_in_ao,ind_p_ao,ind_ga_ao,ind_lig_ao,ind_lig2_ao = get_ao_ind([ind_In,ind_P,ind_Ga,ind_lig,ind_lig2],atoms,nbas,orb_per_atom)

#PYSCF BEGINS	
mol = gto.M()
mol.atom=xyz
mol.basis="def2-svp" #we don't mess around with basis sets here
mol.ecp="def2-svp"
mol.verbose=4
mol.charge=charge
mol.build()

coeff_file_expand=np.load(coeff_file)
mo_mat=coeff_file_expand['arr_0'] # normalized
mo_e = coeff_file_expand['arr_1']

lower_bound=0
upper_bound=nbas
if loc_max > loc_min:
	lower_bound=homo+loc_min+1
	upper_bound=homo+loc_max+1
elif loc_max < loc_min:
	lower_bound=homo-loc_min
	upper_bound=homo-loc_max

if loc_type=="pm" or loc_type=="Pm" or loc_type=="PM" or loc_type=="Pipek" or loc_type=="Pipek Mezey" or loc_type=="Pipek-Mezey":
	mo = Pipek(mol, mo_mat[:,lower_bound:upper_bound+1])
	mo.kernel()
elif loc_type=="boys" or "Boys":
	mo = Boys(mol, mo_mat[:,lower_bound:upper_bound+1])
	mo.kernel()
else:
	print("This type of orbital localization is not yet supported")

for i in range((upper_bound+1)-lower_bound):
	#print(i)
	cubegen.orbital(mol,"mo_"+str(lower_bound+i+1)+"_"+loc_type+".cube",mo.mo_coeff[:,i]) #ideal iso for these is 0.03 to match QChem's 0.02

	#dos analysis
	start=0
	alphas=[]
	for j,atom in enumerate(atoms):
		ao_inds=np.zeros(nbas,dtype=bool)
		norbs=orb_per_atom[atom]
		end=start+norbs
		ao_inds[start:end]=True
		partition=mo.mo_coeff[ao_inds,i]
		alpha=np.sum(np.power([partition],2))
		alphas.append(alpha)
		start=end
	alpha_rank=np.argsort(alphas)
	print("The five major contributions to localized MO", str(lower_bound+i+1), "are", atoms[alpha_rank[-1]], get_atom_specific_index(atoms,alpha_rank[-1]), "with", round(100*alphas[alpha_rank[-1]],2),",",atoms[alpha_rank[-2]], get_atom_specific_index(atoms,alpha_rank[-2]), "with", round(100*alphas[alpha_rank[-2]],2),",",atoms[alpha_rank[-3]], get_atom_specific_index(atoms,alpha_rank[-3]), "with", round(100*alphas[alpha_rank[-3]],2),",",atoms[alpha_rank[-4]], get_atom_specific_index(atoms,alpha_rank[-4]), "with", round(100*alphas[alpha_rank[-4]],2),",",atoms[alpha_rank[-5]], get_atom_specific_index(atoms,alpha_rank[-5]), "with", round(100*alphas[alpha_rank[-5]],2))


