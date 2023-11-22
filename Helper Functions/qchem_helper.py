import numpy as np
import sys

'''
Helper functions for qchem file analysis

python3 compatible

TO DO: combine all into some kind of class or function so that it doesn't have to loop multiple times for each function
'''

def read_lowdin(outfile,atom,index,orb_num,orb_type):
    '''
    Function to extract the MO index corresponding to the first core orbital of a given type on a given atom
    For use in running ddft calculations, mainly

    Inputs: outfile -- name of qchem output file with Lowdin populations (str)
            atom -- atomic species targeted (str)
            index -- atom-specific index of targeted species (str). If 0, just looks for the first atom of the chosen element of any index 
            orb_num -- orbital shell to target. 2 for 2p, for example (int)
            orb_type -- orbital type to target. p for 2p, for example (str). Only works for s,p,d

    Outputs: mo_index -- mo index of target, not 0-indexed (int)
    '''

    if orb_type=="p":
        orb_num=orb_num-1
    elif orb_type=="d":
        orb_num=orb_num-2

    with open(outfile,'r') as out:
        flag=False
        lowdin=[]
        for i,line in enumerate(out):
            if flag==False:
                if line.find("Partial Lowdin")!=-1:
                    flag=True
                    i_start=i
            else:
                if line.find("-----------------------------------------------------------------")!=-1:
                    flag=False
                elif line.find("Populations for Occupied Beta Orbitals")!=-1:
                    flag=False
                    lowdin.pop(-1)
                elif i>i_start+1:
                    line=line[:10]+" "+line[10:]  #gotta add a space to separate "InX"
                    lowdin.append(line.strip().split())

    lowdin.pop(-1)

    orb_count=0
    converged=False
    six_count=0
    orbs=[]
    for i,line in enumerate(lowdin):
        if len(line)==6:
            six_count=six_count+1
        if six_count==1:
            orbs.append(int(line[0]))
        if line[1]==atom and (line[2]==index or index=="0") and line[3]==orb_type:
            for j,pop in enumerate(line[4:]):
                if float(pop)>0.9: #this way it works for the unrestricted girlies too
                    orb_count=orb_count+1
                    if orb_count==orb_num:
                        final_i=i
                        final_j=j
                        converged=True
                        break
        if converged==True:
            break
    n_orbs=max(orbs)

    for i,line in enumerate(lowdin):
        if len(line)==6:
            if final_i-i<n_orbs:
                return line[final_j]



def get_geom_cube(inputfile):
    '''
    Function to extract geometry from qchem cube file

    Inputs:  inputfile  -- name of qchem cube file with geometry
    Outputs: atom_names -- np array with the atom names (str)
             coords     -- np array with the xyz coords, in Bohr (float)
             cube_origin -- np array with the cube origin, which is in Bohr (float)
    '''
    with open(inputfile,'r') as inp:
        line_count = 0
        geom = []
        for i,line in enumerate(inp):
            if line_count==2:
                cube_origin=line.strip().split()
            if line_count>5 and len(line.strip().split()) > 5:
                break
            if line_count > 5:
                geom.append(line.strip().split())
            line_count=line_count+1

    geom = np.array(geom)
    atom_nums = geom[:,0]
    #cube files use numbers for atoms. Currently brute forcing this
    atom_names=[]
    for i,num in enumerate(atom_nums):
        if num=="13":
            atom_names.append("Al")
        elif num=="9":
            atom_names.append("F")
        elif num=="49":
            atom_names.append("In")
        elif num=='15':
            atom_names.append("P")
        elif num=='31':
            atom_names.append("Ga")
        elif num=='30':
            atom_names.append("Zn")
        elif num=='17':
            atom_names.append("Cl")
        elif num=='16':
            atom_names.append("S")
        elif num=='34':
            atom_names.append("Se")
        else:
            print("Element not supported for cube reading (yet)")
    atom_names=np.array(atom_names)
    
    coords = geom[:,2:].astype(float)

    cube_origin=np.array(cube_origin)
    cube_origin=cube_origin[1:].astype(float)

    return atom_names, coords, cube_origin

def get_converged_geom(opt_out):
    '''
    Another Ezra original, right off the cuff 
    Extracts the converged final geometry from a qchem optimization output

    Inputs: opt_out -- .out file with optimization
    Outputs: xyz_coords -- np array with the coordinates of all the atoms (float, NAtoms x 3). Indexed the same as atom_names
             atom_names -- np array with the atom names (str, NAtoms)
    '''
    with open(opt_out,'r') as out:
        flag=0
        lines=[]
        for i,line in enumerate(out):
            if line.find("**  OPTIMIZATION CONVERGED  **")!=-1:
                flag=1
            elif flag==1 and line.find("ATOM                X               Y               Z")!=-1:
                flag=2
            elif flag==2 and line.find("Z-matrix Print:")!=-1:
                break
            elif flag==2:
                lines.append(line.strip().split())
        lines.pop(-1)
        lines=np.array(lines)
        atom_names=lines[:,1]
        xyz_coords=lines[:,2:].astype(float)
        return xyz_coords, atom_names


def my_get_rem_sp(inputfile):
    '''
    Function to extract $rem section and spin/charge from qchem input/output

    Inputs:  inputfile  -- (qchem input/output file)
    Outputs: rem -- string with the rem section separated by \n
             spcharge     -- string with charge and spin
    '''
    with open(inputfile, "r") as inp:
        rem_flag = False
        rem = []
        spcharge = 0
        sp_flag = False
        for i,line in enumerate(inp):
            if sp_flag:
                spcharge = line
                sp_flag = False
            if line.find('$molecule') != -1:
                sp_flag = True
            if line.find('$rem') != -1:
                rem_flag = True
            if rem_flag:
                rem.append(line)
            if rem_flag and line.find('$end') != -1:
                break

        return rem, spcharge