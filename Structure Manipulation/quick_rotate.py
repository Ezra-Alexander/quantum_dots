import numpy as np
import sys
import math
from geom_helper import *

#a quick script to rotate a given .xyz file by a given angle along a given axis

xyz=sys.argv[1] #.xyz
out=sys.argv[2] #name of .xyz to be written
axis=sys.argv[3] #x, y, z or some combination thereof (i.e rotation along the vector {1/sqrt(2),-1/sqrt(2),0} is specified by x-y)
angle=float(sys.argv[4]) #angle of rotation, in degrees

#translate the axis string into a vector
vector=[0,0,0]
neg_flag=False
for i,char in enumerate(axis):
	if char=="x":
		if neg_flag:
			vector[0]=-1
			neg_flag=False
		else:
			vector[0]=1
	elif char=="y":
		if neg_flag:
			vector[1]=-1
			neg_flag=False
		else:
			vector[1]=1
	elif char=="z":
		if neg_flag:
			vector[2]=-1
			neg_flag=False
		else:
			vector[2]=1
	elif char=="-":
		neg_flag=True
	else:
		print("Error with axis variable")
vector=vector/np.linalg.norm(vector)

#read xyzs
coords, atoms = read_input_xyz(xyz)

#convert angle to radians
angle=angle*math.pi/180

#rotation matrix
rot_mat=np.array([[ math.cos(angle)+(vector[0]**2*(1-math.cos(angle))) , (vector[0]*vector[1]*(1-math.cos(angle)))-(vector[2]*math.sin(angle)) , (vector[0]*vector[2]*(1-math.cos(angle)))+(vector[1]*math.sin(angle)) ] , [ (vector[0]*vector[1]*(1-math.cos(angle)))+(vector[2]*math.sin(angle)) , math.cos(angle)+(vector[1]**2*(1-math.cos(angle))) , (vector[2]*vector[1]*(1-math.cos(angle)))-(vector[0]*math.sin(angle)) ] , [ (vector[0]*vector[2]*(1-math.cos(angle)))-(vector[1]*math.sin(angle)) , (vector[2]*vector[1]*(1-math.cos(angle)))+(vector[0]*math.sin(angle)) , math.cos(angle)+(vector[2]**2*(1-math.cos(angle))) ]])


#do rotation
new_coords=[]
for i,atom in enumerate(atoms):
	coord_vector=np.transpose(coords[i])
	rotated_coords=np.matmul(rot_mat,coord_vector)
	new_coords.append(rotated_coords)
new_coords=np.array(new_coords)

write_xyz(out,atoms,new_coords)

