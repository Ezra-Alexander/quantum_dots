import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from qd_helper import *
import copy
from geom_helper import *

#just a quick script to get a rms distortion of an atom in a QD before and after a defect

new_xyz=sys.argv[1] #the new xyz file in question
old_xyz=sys.argv[2] #the xyz file before the defect
element=sys.argv[3] #the element of the specified atom
index=int(sys.argv[4]) #element specific index of target. We assume this element-specific index doesn't change

cutoff = 2.8 #cutoff for 2 atoms to be considered bonded. Arbitrary and may need to be changed

new_coords,new_atoms = read_input_xyz(new_xyz)
old_coords,old_atoms = read_input_xyz(old_xyz)

target = 0
overall_index_new="bad"
for i,atom in enumerate(new_atoms):
	if atom == element:
		target = target+1
		if target == index:
			overall_index_new=i

target = 0
overall_index_old=0
for i,atom in enumerate(old_atoms):
	if atom == element:
		target = target+1
		if target == index:
			overall_index_old=i

rmsd=get_rms_distortion(new_coords,old_coords,overall_index_new, overall_index_old,cutoff)

print("Root mean square distortion for", element, index, "is ", round(rmsd,3))