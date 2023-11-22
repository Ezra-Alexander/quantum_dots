# from qd_helper import *
import copy
import numpy as np
import matplotlib.pyplot as plt
import math

import numpy as np
from matplotlib import pyplot as plt

#Including only functions which I have written. So, much of the code may not actually work 

def bring_to_4c(atoms,coords,target_index, ligand,threshold):
    '''
    A function that passivated a given atom from as low as 1c up to four-coordinate in a tetrahedral geometry with a chosen ligand

    Inputs:
        atoms: np array of atom names, in order. NAtoms long
        coords: np array of atomic coords, in the same order. NAtoms x 3
        target_index: int. the overall index, python adjusted, of the atom you are targeting
        ligand: str. the element you will use to passivate the target. Note that bond lengths are hard-coded
        threshold: float. the maximum distance between two atoms for them to be considered bonded
    Outputs:
        temp_atoms: updated np array of atom names, in order. NAtoms long
        temp_coords: updated np array of atomic coords, in the same order. NAtoms x 3
    '''

    bond_lengths={"In":{"F":2.02,"Cl":2.37},"Ga":{"F":1.79,"Cl":2.19},"P":{"O":1.53}} #hard coded (terminal) bond lengths. not everything will be supported. Don't have to be perfect, these are just pre-opt estimates

    bond_length=bond_lengths[atoms[target_index]][ligand]

    temp_coords=np.copy(coords)
    temp_atoms=np.copy(atoms)

    dists=dist_all_points(temp_coords)
    connectivity=connectivity_finder(dists,threshold)
    if len(connectivity[target_index])==1:
        #the hard part for a 1-coordinate target is choosing the direction of the new bond
        anchor=connectivity[target_index][0]
        anchor_bonds=connectivity[anchor]
        vector=temp_coords[anchor]-temp_coords[anchor_bonds[0]]
        unit=vector/np.linalg.norm(vector)
        added_coords=temp_coords[target_index]+(unit*bond_length)
        temp_coords=np.append(temp_coords,[added_coords],axis=0)
        temp_atoms=np.append(temp_atoms,ligand)

    dists=dist_all_points(temp_coords)
    connectivity=connectivity_finder(dists,threshold)
    if len(connectivity[target_index])==2:
        v1=temp_coords[target_index]-temp_coords[connectivity[target_index][0]]
        v2=temp_coords[target_index]-temp_coords[connectivity[target_index][1]]

        parallel=(v1+v2)/np.linalg.norm(v1+v2)
        perp=np.cross(v1,v2)
        perp_unit=perp/np.linalg.norm(perp)

        adjacent=bond_length*math.cos(math.radians(109.5/2))
        opposite=bond_length*math.sin(math.radians(109.5/2))

        added_coords=temp_coords[target_index]+(parallel*adjacent)+(perp_unit*opposite)
        temp_coords=np.append(temp_coords,[added_coords],axis=0)
        temp_atoms=np.append(temp_atoms,ligand)


    dists=dist_all_points(temp_coords)
    connectivity=connectivity_finder(dists,threshold)
    if len(connectivity[target_index])==3:
        temp_atoms,temp_coords=geom_adder(temp_atoms,temp_coords,target_index,ligand)

    return temp_atoms, temp_coords



def geom_adder(atoms,coords,target_index,add_element):
    '''
    From the script of the same name
    A script meant to add a chosen element to a given 3-coordinate target in the 4th tetrahedral position
    Technically also works for adding something to a 4+ coordinate target, but it ends up putting it somewhere weird usually

    Inputs:
        atoms: np array of atom names, in order. NAtoms long
        coords: np array of atomic coords, in the same order. NAtoms x 3
        target_index: int. the overall index, python adjusted, of the atom you are targeting
        add_element: str. the element you wish to add. note that target-add bond lengths are hard coded and not all may be supported yet
    '''
    #find vector along 4th dimension
    ind_In = atoms=='In'
    ind_P = atoms=='P'
    ind_Ga = atoms=="Ga"
    ind_InGa = np.logical_or(ind_In,ind_Ga)
    ind_cat=np.logical_or(ind_InGa,atoms=="Al")
    ind_InP= np.logical_or(ind_In,ind_P)
    ind_lig = (atoms == "Cl")

    all_dists,inp_dists,inf_dists,inpf_dists,pin_dists = get_dists(coords,ind_cat,ind_P,ind_lig)
    #all_dists,inp_dists,inf_dists,inpf_dists,pin_dists = get_dists(coords,ind_Ga,ind_P,ind_lig)

    target_dists = all_dists[target_index]
    #sorted_dists = np.sort(target_dists)
    indexes = [1000,1000,1000]
    dists = [1000,1000,1000]
    for i,dist in enumerate(target_dists):
        if dist < max(dists) and i!=target_index:
            indexes.pop(dists.index(max(dists)))
            dists.pop(dists.index(max(dists)))
            indexes.append(i)
            dists.append(dist)

    ind_1 = indexes[0]
    ind_2 = indexes[1]
    ind_3 = indexes[2]

    #bond_1 = sorted_dists[1]
    #bond_2 = sorted_dists[2]
    #bond_3 = sorted_dists[3]

    #this assumes that the three closest dists are not exactly equal
    #ind_1 = np.where(target_dists==bond_1)[0][0]
    #ind_2 = np.where(target_dists==bond_2)[0][0]
    #ind_3 = np.where(target_dists==bond_3)[0][0]

    # dummy = (coords[ind_1] + coords[ind_2] + coords[ind_3])/3

    #dummy = (coords[ind_2] + coords[ind_3])/3 #temp

    v1=coords[target_index]-coords[ind_1]
    v2=coords[target_index]-coords[ind_2]
    v3=coords[target_index]-coords[ind_3]

    u1=v1/np.linalg.norm(v1)
    u2=v2/np.linalg.norm(v2)
    u3=v3/np.linalg.norm(v3)

    vector=u1+u2+u3

    #vector = coords[target_index]-dummy
    mag = np.linalg.norm(vector)
    unit = vector/mag

    #Need to add bond lengths of each type as I go
    if add_element == "In" and atoms[target_index] == "P":
        new_coords = coords[target_index] + unit*2.58
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element == "F" and atoms[target_index] == "In":
        new_coords = coords[target_index] + unit*2.2
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element == "O" and atoms[target_index] == "In":
        new_coords = coords[target_index] + unit*2.05
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element == "Ga" and atoms[target_index] == "P":
        new_coords = coords[target_index] + unit*2.36
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="F" and atoms[target_index]=="Ga":
        new_coords = coords[target_index] + unit*2.12
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="Al" and atoms[target_index]=="P":
        new_coords = coords[target_index] + unit*2.32
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="Cl" and atoms[target_index]=="In":
        new_coords = coords[target_index] + unit*2.36
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="Cl" and atoms[target_index]=="Zn":
        new_coords = coords[target_index] + unit*2.186
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="Cl" and atoms[target_index]=="Al":
        new_coords = coords[target_index] + unit*2.116
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="Cl" and atoms[target_index]=="Ga":
        new_coords = coords[target_index] + unit*2.19
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    elif add_element=="O" and atoms[target_index]=="P":
        new_coords = coords[target_index] + unit*1.53
        print("Done!")
        new_atoms = np.append(atoms,add_element)
        final_coords = np.append(coords,[new_coords], axis=0)
    else:
        print("Warning! The bond length between the element you are adding and the element you are targeting is not yet supported!")
    return new_atoms,final_coords

def get_label(atoms,coords,target_element,target_index,threshold,ligand):
    '''
    Function that uses the connectivity mapping to label a given atom as being either under-coordinated, surface, or bulk
    With delta-DFT in mind but there are probably other uses

    Inputs:
        atoms: np array of atom names, in order. NAtoms long
        coords: np array of atomic coords, in the same order. NAtoms x 3
        target_element: string. the element of the atom you are targeting. this could be written with overall indices but i prefer this
        target_index: int. the atom-specific index of the atom you are targeting
        threshold: float. the upper bound for interatomic distance for two atoms to be considered "bonded"
    Outputs:
        label: string. either "Under-Coordinated", "Bulk", "Bound to Trap", or "Bound to Dopant"
            Anything other than In, Ga, Al, P, or ligand is considered a dopant
            Only bulk "near dopants" will be included
    '''
    dists=dist_all_points(coords)
    connectivity=connectivity_finder(dists,threshold)
    count=0
    uc_flag=False
    dopant_flag=False
    for i,atom in enumerate(atoms):
        if atom==target_element:
            count=count+1
            if count==target_index:
                if len(connectivity[i])<4:
                    # for j,bond in enumerate(connectivity[i]):
                        # if len(connectivity[bond])<4 and atoms[bond]!=ligand:
                        #     print("This shouldn't be happening")
                        #     return "Bound to Trap"
                    return "Under-Coordinated"
                else:
                    for j,bond in enumerate(connectivity[i]):
                        if len(connectivity[bond])<4 and atoms[bond]!=ligand:
                            uc_flag=True
                            #return "Bound to Trap"
                        if atoms[bond]!=ligand and atoms[bond]!="In" and atoms[bond]!="Ga" and atoms[bond]!="Al" and atoms[bond]!="P":
                            dopant_flag=True
                            #return "Bound to Dopant"
                        # elif atoms[bond]==ligand:
                        #     return "Surface"
                        # else:
                        #     for k,bond2 in enumerate(connectivity[bond]):
                        #         if atoms[bond2]==ligand:
                        #             return "Surface"
                    if dopant_flag==True:
                        return "Bound to Dopant"
                    elif uc_flag==True:
                        return "Bound to Trap"


                    return "Bulk"



def to_atom_specific_index(atoms,overall_index):
    '''
    Function that converts from overall indices (i.e. python index 239 is the 240th atom in the structure)
    to atom-specific indices (i.e. P 12, the 12th Phosphorus)

    Inputs:
        atoms: np array of atom names, in order. shape (Natoms)
        overall_index: positive integer, 0-based indexed.
    '''
    target_atom=atoms[overall_index]
    count=0
    for i,atom in enumerate(atoms):
        if atom==target_atom:
            count=count+1
            if i==overall_index:
                return count


def get_underc_index_variable(xyz,ind_Cd,ind_Se,ind_lig,lig_cutoff,p_cutoff,nncutoff):
    '''
    Function that finds undercoordinated Cd and Se atoms in a QD.

    Inputs:
        xyz: np array of xyz coordinates for the QD. shape (Natoms,3)
        ind_Cd: boolean array of shape Natoms, indexing the Cd atoms
        ind_Se: boolean array of shape Natoms, indexing the Se atoms
        ind_lig: boolean array of shape Natoms, indexing the ligand atoms
        lig_cutoff: upper cutoff for In-F bonds
        p_cutoff: upper cutoff for In-P bonds
        nncutoff: number of nearest neighbors to be considered "fully coordinated"
                  (< this classified as "undercoordinated")
    '''
    all_dists,cdse_dists,cdlig_dists,cdselig_dists,secd_dists = get_dists(xyz,ind_Cd,ind_Se,ind_lig)

    all_nn=num_nn_variable(ind_Cd,ind_Se,ind_lig,all_dists,lig_cutoff,p_cutoff)

    cd_underc_ind=[]
    se_underc_ind=[]
    for i,n_nn in enumerate(all_nn):
        if ind_Cd[i]:
            cd_underc_ind.append(n_nn<nncutoff)
        elif ind_Se[i]:
            se_underc_ind.append(n_nn<nncutoff)

    cd_underc_ind=np.array(cd_underc_ind)
    se_underc_ind=np.array(se_underc_ind)

    return cd_underc_ind,se_underc_ind


def num_nn_variable(ind_Cd,ind_Se,ind_lig,dist_list,lig_cutoff,p_cutoff):
    '''
    Honestly not sure if Lexie already implemented this or not

    Does the coordination # counting, but uses different cutoffs for cation-ligand bonds and cation-anion bonds

    Inputs:
        ind_Cd: np array of which atoms are our cation, Natoms long
        ind_Se: np array of which atoms are our anion, Natoms long
        ind_lig: np array of which atoms are our ligand, Natoms long
        dist_list: Natoms x Natoms, the distance of each atom to each other atom. Includes self-distances
        lig_cutoff: the upper bond length cutoff for cation-ligand bonds, Angstroms
        p_cutoff: the upper bond length cutoff for anion-ligand bonds, Angstroms
    
    Outputs:
        nn_list: an array of the number of nearest neighbors for each atom. Size (Natoms,).
                 nn_list[i] = # of nearest neighbors for atom i
    '''
    nn_list=[]
    for i,atom in enumerate(dist_list):
        nn_count=0
        for j,dist in enumerate(atom):
            if i!=j: #this is how we're dealing with the self-bonds
                if ind_Cd[i]:
                    if ind_lig[j]:
                        if dist<lig_cutoff:
                            nn_count=nn_count+1
                    else:
                        if dist<p_cutoff:
                            nn_count=nn_count+1
                elif ind_Se[i]:
                    if dist<p_cutoff:
                            nn_count=nn_count+1
                elif ind_lig[i]:
                    if dist<lig_cutoff:
                            nn_count=nn_count+1
                else: #going to use the bigger cutoff for anything else
                    if dist<p_cutoff:
                            nn_count=nn_count+1
        nn_list.append(nn_count)

    nn_list=np.array(nn_list)


    return nn_list

def get_rms_distortion(new_coords,old_coords,target_index_new,target_index_old,cutoff):
    '''
    The goal of this function is to compute a sort of root-mean-squared deviation between the structure of an atom before and after a certain defect

    Centers the target atom in the same place and computes the difference between its bonded constituents

    Obviously won't work if the atom changes in coordination #
    Too lazy to code this to work with a change in the atom's indices

    Inputs:
        new_coords: np array of the coordinates of all atoms in the final structure, the one we are interested in
        old_coords: np array of the original (but still optimized) coordinates of all atoms before said defect
        target_index_new: integer, the overall index of the atom we are interested in in the new xyz
        target_index_old: integer, the overall index of the atom we are interested in in the old xyz
        cutoff: integer, a manual parameter that specifies what counts as "bonded"
    Outputs:
        rmsd: float, the rms deviation of the positions of all bonded atoms to our target
    '''

    new_dists=dist_all_points(new_coords) #all distances in new dot
    old_dists=dist_all_points(old_coords) #all distances in old dot

    #recenter coords
    new_cent_coords=new_coords-new_coords[target_index_new]
    old_cent_coords=old_coords-old_coords[target_index_old]

    target_new_dists=new_dists[target_index_new]
    target_old_dists=old_dists[target_index_old]

    bonded_indeces_new=[]
    for i,dist in enumerate(target_new_dists):
        if dist < cutoff and dist > 0:
            bonded_indeces_new.append(i)

    bonded_indeces_old=[]
    for i,dist in enumerate(target_old_dists):
        if dist < cutoff and dist > 0:
            bonded_indeces_old.append(i)


    #we now need to "match" the bonded atoms with themselves. We also get the distances for free
    deltas=[]
    for i,bond in enumerate(bonded_indeces_new):
        diffs=[]
        for j,old_bond in enumerate(bonded_indeces_old):
            diffs.append(np.sqrt(np.sum((new_cent_coords[bond] - old_cent_coords[old_bond])**2)))
        deltas.append((min(diffs)))
    deltas=np.array(deltas)

    #now compute the RMS
    num=np.sum(deltas**2)
    rms=np.sqrt(num/4)
    return rms



def connectivity_finder(dists,cutoff):
    '''
    The people have spoken and they say "more ezra code"

    makes an NAtoms list containing, for each atom, the indices of all other atoms within cutoff of that atom

    Inputs:
        dists: np array with distances between all atoms, size (Natoms, Natoms)
        cutoff: cutoff for a nearest neighbor distance
    Outputs:
        connect: NAtom alist of NBonded lists of python-adjusted indices for all other atoms within cutoff of that atom
    '''
    connectivity=[]
    for i,atom in enumerate(dists):
        bonds=[]
        for j,dist in enumerate(atom):
            
            if float(dist) < cutoff and i!=j:
                bonds.append(j)
        connectivity.append(bonds)

    return connectivity





def get_coplane(coords,dists,target_index,cutoff):
    '''
    Another iconic piece of ezra code

    finds a "coplanarity metric", i.e. some sort of distance from a plane (or line) for a given atom in a structure

    Inputs:
        coords: np array of xyz coordinates for the QD. Has shape (Natoms,3)
        dists: np array with distances between all atoms, size (Natoms, Natoms)
        target_index: the overall index of the atom in the .xyz file, adjusted to base 0 for python
        cutoff: cutoff for a nearest neighbor distance
    Outputs:
        the coplanarity metric. Defined differently for 2c, 3c, and 4c atoms
    '''
    bonded_indeces = []
    for i,atom in enumerate(coords):
        if dists[target_index][i] < cutoff and i != target_index: # and atom!="F": #for computing In/Ga-4c
            bonded_indeces.append(i)

    if len(bonded_indeces) > 3:
        #if len(bonded_indeces) == 5:
        #    print(target_index, "is 5 CMC")
        #print("Hope this is supposed to be 4-coordinate:")

        umbrella=1
        min_coplane = 10 #arbitrary high number
        for i,bond in enumerate(bonded_indeces): #loop through bonded atoms, find coplane metric excluding that atom
            v1 = coords[bonded_indeces[i-1]] - coords[bonded_indeces[i-2]]
            v2 = coords[bonded_indeces[i-3]] - coords[bonded_indeces[i-2]]
            norm_v = np.cross(v1,v2) #plane equation comes from normal vector to plane
            a = norm_v[0]
            b = norm_v[1]
            c = norm_v[2]
            d = -a*coords[bonded_indeces[i-1]][0]-b*coords[bonded_indeces[i-1]][1]-c*coords[bonded_indeces[i-1]][2]
            dist_from_plane = (a*coords[target_index][0]+b*coords[target_index][1]+c*coords[target_index][2]+d)/math.sqrt(a**2 + b**2 + c**2)

            if dist_from_plane<0: #sometimes it'll be negative
                dist_from_plane = dist_from_plane*-1

            #print(dist_from_plane, bond)

            if dist_from_plane < min_coplane:
                min_coplane=dist_from_plane

                i_dist=np.linalg.norm(coords[bonded_indeces[i]]-coords[target_index]) #dist of 4th point from center
                i_dist_from_plane=abs((a*coords[bonded_indeces[i]][0]+b*coords[bonded_indeces[i]][1]+c*coords[bonded_indeces[i]][2]+d)/math.sqrt(a**2 + b**2 + c**2))
                other_d1=np.linalg.norm(coords[bonded_indeces[i]]-coords[bonded_indeces[i-1]])
                other_d2=np.linalg.norm(coords[bonded_indeces[i]]-coords[bonded_indeces[i-2]])
                other_d3=np.linalg.norm(coords[bonded_indeces[i]]-coords[bonded_indeces[i-3]])
                if (i_dist > i_dist_from_plane) and (max(other_d1,other_d2,other_d3)<1.5*i_dist): #trying to distingish structures with umbrella like arrangements from the normal tetrahedron
                    umbrella=-1
                    #print(max(other_d1,other_d2,other_d3))

        return min_coplane*umbrella

    elif len(bonded_indeces) < 3:
        print("Hope this is supposed to be 2-coordinate:")

        v1 = coords[bonded_indeces[1]] - coords[bonded_indeces[0]]
        u1 = v1/np.linalg.norm(v1)
        d1=coords[target_index]-coords[bonded_indeces[0]]
        dist_vect=d1-np.dot(d1,u1)*u1
        dist_from_line=np.linalg.norm(dist_vect)

        return dist_from_line

    else:
        v1 = coords[bonded_indeces[1]] - coords[bonded_indeces[0]]
        v2 = coords[bonded_indeces[2]] - coords[bonded_indeces[0]]
        norm_v = np.cross(v1,v2) #plane equation comes from normal vector to plane
        a = norm_v[0]
        b = norm_v[1]
        c = norm_v[2]
        d = -a*coords[bonded_indeces[1]][0]-b*coords[bonded_indeces[1]][1]-c*coords[bonded_indeces[1]][2]
        dist_from_plane = (a*coords[target_index][0]+b*coords[target_index][1]+c*coords[target_index][2]+d)/math.sqrt(a**2 + b**2 + c**2)

        if dist_from_plane<0:
            dist_from_plane = dist_from_plane*-1

        return dist_from_plane


def get_off_tet_index(QD_xyz, ind_In, ind_P,ind_lig,in_underc_ind,p_underc_ind,dists,cutoff,angle_cutoff):
    '''
    This one's fully written by me

    Function that finds strongly off-tetrahedral 4 (and 5) Coordinate In and P atoms in a QD
    Defined as any 4-c P or In with at least three bond angles centered at that atom above angle_cutoff away from 109.5 degrees

    Inputs:
        xyz: np array of xyz coordinates for the QD. Has shape (Natoms,3)
        ind_In: boolean array of shape Natoms, indexing the In atoms
        ind_P: boolean array of shape Natoms, indexing the P atoms
        ind_lig: boolean array of shape Natoms, indexing the ligand atoms
        in_underc_ind: boolean array of shape NIndium, indexing the < 4-c In atoms
        p_underc_ind: boolean array of shape NPhosphorous, indexing the < 4-c P atoms 
        dists: np array with distances between all atoms, size (Natoms, Natoms)
        cutoff: cutoff for a nearest neighbor distance
        angle_cutoff: minimum degrees from 109.5 for an angle to be considered "off-tetrahedral"
    Outputs:
        in_off_tet_ind: boolean array of shape NIndium, indexing off-tetrahedral In-4c
        p_off_tet_ind: boolean array of shape NPhosphorous, indexing off-tetrahedral P-4c
        in_angles: list with angles centered on each In. Size (NIndium, Nangles). Unchanged from get_angles
        p_angles: list with angles centered on each P. Size (NPhosphorous, Nangles). Unchanged from get_angles
    '''

    in_angles,p_angles = get_angles(QD_xyz,ind_In, ind_P,dists,cutoff)

    in_off_tet_ind = []
    for i,indium in enumerate(in_angles):
        ok = True
        count = 0
        if in_underc_ind[i]:
            ok = False
            in_off_tet_ind.append(False)
        for j,angle in enumerate(indium):
            if abs(angle-109.5)>=angle_cutoff and ok and count < 1:
                count = count+1
            elif abs(angle-109.5)>=angle_cutoff and ok and count == 1:
                in_off_tet_ind.append(True)
                ok = False
            if j+1 == len(indium) and ok:
                in_off_tet_ind.append(False)                

    p_off_tet_ind = []
    for i,phos in enumerate(p_angles):
        ok = True
        count = 0
        if p_underc_ind[i]:
            ok = False
            p_off_tet_ind.append(False)
        for j,angle in enumerate(phos):
            if abs(angle-109.5)>=angle_cutoff and ok and count < 1:
                count = count+1
            elif abs(angle-109.5)>=angle_cutoff and ok and count == 1:
                p_off_tet_ind.append(True)
                ok = False
            if j+1 == len(phos) and ok:
                p_off_tet_ind.append(False) 

    #print(p_off_tet_ind)


    return in_off_tet_ind, p_off_tet_ind, in_angles, p_angles

def get_angles(QD_xyz,ind_In,ind_P,dists,cutoff):
    '''
    Another ezra classic

    Function that calculates all angles centered at each In and P atom between atoms within cutoff

    Inputs:
        QD_xyz: np array of xyz coordinates for the QD. Has shape (Natoms,3)
        ind_In: boolean array of shape Natoms, indexing the In atoms
        ind_P: boolean array of shape Natoms, indexing the P atoms
        ind_lig: boolean array of shape Natoms, indexing the ligand atoms
        dists: np array with distances between all atoms, size (Natoms, Natoms)
        cutoff: cutoff for a nearest neighbor distance
    Outputs: 
        in_angles: list with angles centered on each In. Size (NIndium, Nangles for that atom (depends on number of nearest neighbors))
        p_angles: list with angles centered on each P. Size (NPhosphorous, Nangles)
    '''

    #get all angles
    angles=[]
    for i,center in enumerate(QD_xyz):
        atom = []
        for j,atom1 in enumerate(QD_xyz):
            if dists[i][j]<cutoff and i!=j:
                for k,atom2 in enumerate(QD_xyz):
                    if dists[i][k]<cutoff and i!=k and j!=k and k>j:
                        dist1 = dists[i][j]
                        dist2 = dists[i][k]
                        dist3 = dists[j][k]
                        num = dist1**2 + dist2**2 - dist3**2
                        denom = 2*dist1*dist2
                        angle = math.acos(num/denom)
                        atom.append(math.degrees(angle))
        angles.append(atom)
    #all_angles = np.array(angles)

    in_angles = []
    p_angles = []
    for i,angle in enumerate(angles):
        if ind_In[i]:
            in_angles.append(angle)
        if ind_P[i]:
            p_angles.append(angle)


    return in_angles, p_angles

def get_off_pairs(dists,p_off_tet_ind,pp_cutoff,ind_P):
    '''
    written by ez

    Function that finds "adjacent" pairs of off-tetrahedral P (within pp_cutoff)

    Inputs:
        dists: np array with distances between all atoms, size (Natoms, Natoms)
        p_off_tet_ind: boolean array of shape NPhosphorous, indexing off-tetrahedral P-4c
        pp_cutoff: cutoff for two P to be considered adjacent, in Angstroms
        ind_P: boolean array of shape Natoms, indexing the P atoms
    Outputs:
        p_pairs:
    '''

    p_pairs = []
    p_count= -1
    #print(len(p_off_tet_ind))
    for i,atom in enumerate(dists):
        if ind_P[i]:
            #print("i",i)
            p_count=p_count+1
            if p_off_tet_ind[p_count]:
                p_count2 = -1
                for j,atom2 in enumerate(atom):
                    if ind_P[j]:
                        #print("j",j)
                        p_count2=p_count2+1
                        if p_off_tet_ind[p_count2] and atom2 < pp_cutoff and i != j and j > i:
                            p_pairs.append([p_count+1,p_count2+1])

    return p_pairs

