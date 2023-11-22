import numpy as np
from matplotlib import pyplot as plt

def get_atom_specific_index(atoms,overall_index):
    '''
    Yet another ezra classic
    Take a list of atoms (from an .xyz file) and an overall index (i.e. the 70th atom overall) and return the atom specific index (i.e. P 25)
    Inputs: atoms -- np array or list of atom names (str)
            overall_index -- integer overall atomic index, already 0-based
    Outputs: atom_specific_index -- integer atom specific index, 1-based for human comprehension
    '''
    target_atom=atoms[overall_index]
    target_count=0
    for i,atom in enumerate(atoms):
        if atom==target_atom:
            target_count=target_count+1
            if i==overall_index:
                return target_count