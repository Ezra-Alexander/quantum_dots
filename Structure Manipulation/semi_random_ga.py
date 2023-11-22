import sys
from qchem_helper import read_xyz,write_input,write_xyz
import numpy as np

#for making semi-random InGaP dots from an InP dot. It can be used to sub other things for In/Ga into InGaP as well 
#2 modes. First is to substitute the ones specified
#specify the In indices of any number of In you want to substitute
#if under half the In, the rest of the substitutions are done randomly
#Second mode is substitute anything but the ones specified

#THIS CODE ASSUMES A CHARGE OF 0 (not that it would be an impossible fix)

xyz = sys.argv[1] #the one you're taking in
output = sys.argv[2] #the name of the opt_qchem.in file you're writing
sub = sys.argv[3] #should either be "sub" or "nosub" for modes 1 or 2 respectively
dopant=sys.argv[4] #in case you want something other than Ga substituted
num_ga=int(sys.argv[5]) #number of requested substitutions 

indeces=sys.argv[6:] #the In indices of the In you want guaranteed substituted / not substituted
indeces=[int(x)-1 for x in indeces] 
indeces.sort(reverse=True)

#read .xyz
coords, atoms = read_xyz(xyz)
natoms = len(atoms)
#in_ind = atoms=='In'
in_ind = np.logical_or(atoms=='In',atoms=='Ga') # a temporary measure, when substituting not Ga into InGaP
num_in = np.count_nonzero(in_ind)


if sub == "sub":
	#do non-random substitutions
	in_to_ga = np.where(in_ind)[0][indeces]
	atoms[in_to_ga] = dopant

	#re-count In
	#in_ind = atoms=='In'
	in_ind = np.logical_or(atoms=='In',atoms=='Ga') # a temporary measure
	num_in = np.count_nonzero(in_ind)

	#disclaimer
	if len(indeces) > num_ga:
		print("Over half of the In will be substituted")
	elif len(indeces) == num_ga:
		print("No random substitutions performed - exactly half specified")
	else: #do the random subs
		num_new_ga = num_ga - len(indeces)
		random_in = np.random.choice(num_in,num_new_ga,replace=False)
		in_to_ga = np.where(in_ind)[0][random_in]
		atoms[in_to_ga] = dopant
elif sub=="nosub":
	if len(indeces) > (num_in/2)+0.5:
		print("You've requested that over half of the Indium be frozen")

	indium=np.arange(num_in)
	for i,ind in enumerate(indeces):
		if ind in indium:
			indium=np.delete(indium,np.where(indium==ind))

	random_in = np.random.choice(indium,num_ga,replace=False)
	in_to_ga = np.where(in_ind)[0][random_in]
	atoms[in_to_ga] = dopant

else:
	print("3rd variable should be either 'sub' to guarantee substitution of the specified indices, or 'nosub' to guarantee the specified indeces are not substituted")
	sys.exit()


rem = '$rem\njobtype opt\nmethod pbe\nbasis def2-svp\necp def2-ecp\nmem_total 5000\nmem_static 1000\nithrsh_dft 15\nscf_algorithm rca\nscf_convergence 6\nmax_scf_cycles 500\ngeom_opt_update 5 \n $end'
sp = '0 1\n'
#write_input(output,atoms,coords,rem,sp)
write_xyz(output,atoms,coords)
