#This code is meant to take all but certain indices of atoms and replace them with a different atom (Nitrogen) so i can make nice pictures in vesta
import numpy as np
import sys


def read_xyz(input_xyz):
    '''
    Function that reads xyz file into arrays

    Inputs: input_xyz -- .xyz file with the QD coordinates
    Outputs: xyz_coords -- np array with the coordinates of all the atoms (float)
             atom_names -- np array with the atom names (str)
    '''
    xyz_coords = np.loadtxt(input_xyz,skiprows=2,usecols=(1,2,3))
    atom_names = np.loadtxt(input_xyz,skiprows=2,usecols=(0,),dtype=str)
    return xyz_coords, atom_names

def write_geom(xyz_file,atom_names,atom_xyz):
    '''
    Function to write the geometry part of an xyz/input file

    Inputs: xyz_file   -- file object, open for writing, to write xyz to
            atom_names -- np array or list of atom names (str)
            atom_xyz   -- np array or list of atom xyz coordinates (float)

    Outputs: Writes to out_file
    '''

    for i, atom in enumerate(atom_names):
        xyz_file.write('{:2}     {:15.10f}     {:15.10f}    {:15.10f}\n'.format(atom, atom_xyz[i][0], atom_xyz[i][1], atom_xyz[i][2]))
    return

def write_xyz(out_file, atom_names, atom_xyz, comment=''):
    '''
    Function that writes xyz coordinate arrays to .xyz file

    Inputs: out_file   -- name of the file to write coordinates to
            atom_names -- np array or list of atom names (str)
            atom_xyz   -- np array or list of atom xyz coordinates (float)
            comment    -- comment line to write in the xyz file

    Outputs: Writes to xyz_file
    '''
    with open(out_file,'w') as xyz_file:
        xyz_file.write(str(len(atom_names))+'\n')
        xyz_file.write(comment+'\n')
        write_geom(xyz_file,atom_names,atom_xyz)
    return


xyz_file = sys.argv[1]

xyz,atoms=read_xyz(xyz_file)

ind = []
phos = []
oxy = []
carb = []
hyd = []

for i,atom in enumerate(atoms):
	if atom == "In":
		ind.append([atom,i])
	if atom == "P":
		phos.append([atom,i])
	if atom == "O":
		oxy.append([atom,i])
	if atom == "C":
		carb.append([atom,i])
	if atom == "H":
		hyd.append([atom,i])


#this section needs to be edited for each run
for i,atom in enumerate(oxy):
	if i != 55 and i != 97 and i != 98 and i != 70 and i != 99 and i != 53:
		atom[0] = "N"

for i,atom in enumerate(carb):
	if i != 54 and i != 57 and i != 67 and i != 65 and i != 50 and i != 44:
		atom[0] = "N"


for i,atom in enumerate(hyd):
	if i != 151 and i != 150 and i != 27 and i != 110 and i != 108 and i != 107 and i != 112 and i!=34 and i != 94 and i != 26 and i != 97:
		atom[0] = "N"

changed_atoms = ind+phos+oxy+carb+hyd

changed_atoms = sorted(changed_atoms, key = lambda x: x[1])


out_atoms = []

for atom in changed_atoms:
	out_atoms.append(atom[0])

write_xyz("picture_perfect.xyz", out_atoms, xyz)

