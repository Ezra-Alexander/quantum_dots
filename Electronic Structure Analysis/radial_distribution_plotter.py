import numpy as np
import numpy.linalg as npl
import sys
from matplotlib import pyplot as plt
import cube_tools
from qchem_helper import *

# a (hopefully) simple script to compute and plot a radial probability distribution for one or more MOs in a QD from cube files

#xyz=sys.argv[1] #the qd structure. must be careful that it matches the cube file
cube_files=indices=[x for x in sys.argv[1:]] #any number of cube files

dr=40/150 #comes from the paramters for my cube files. 40 angstroms, 150 grid points in each dimension

#compute center of mass of the dot
#info the code has
masses = {"Ga":69.723, "P":30.974, "In":114.82,"Cl":35.453,"F":18.998403162,"Al":26.982,"Zn":65.38,"Se":78.97,"S":32.06}

#if xyz[-4:]==".xyz":
#	coords,atoms=read_xyz(xyz)
#elif  xyz[-3:]==".in":
#	atoms,coords=get_geom_io(xyz)
#	rem,sp=my_get_rem_sp(xyz)
#	new_rem=""
#	for i,line in enumerate(rem):
#		new_rem=new_rem+line
#else:
#	print("File type not supported")

atoms,coords,cube_origin=get_geom_cube(cube_files[0])

#choose center of dot
#compute center of mass
cent =np.array([0.0])
mass = 0
for i,atom in enumerate(atoms):
	cent = cent + coords[i]*masses[atom]
	mass = mass + masses[atom]
	com = cent/mass
#print(com)
#print(cube_origin)
com=com-cube_origin
com=com*0.529177210903 #because qchem cube files are written in Bohr, not Angstroms

#read cubes
cubes=[]
for i,cube in enumerate(cube_files):
		cubes.append(cube_tools.cube(cube))

#do the integration
#radii_grid = np.arange(0,30,dr)
radii_grid, dr = np.linspace(0., 20, 75, retstep=True, endpoint=False)
radial_distribution = np.zeros(len(radii_grid))
for i,radius in enumerate(radii_grid):
	for j,cube in enumerate(cube_files):
		rCube=cubes[j]
		#print(rCube.cube_int_ref(com, radius+dr))
		radial_int = ((rCube.cube_int_ref(com, radius+dr) - rCube.cube_int_ref(com, radius))/dr)**2
		radial_distribution[i]=radial_distribution[i]+radial_int
	#if radius>0:
	#	radial_distribution[i]=radial_distribution[i]/(radius**2) 



#plot it
plt.figure()
plt.plot(radii_grid,radial_distribution,'C0')
plt.xlim(0,20)
plt.ylim(0,max(radial_distribution)*1.2)
plt.stem([8,11],[max(radial_distribution)*1.2,max(radial_distribution)*1.2],'Black') # put a stem at the interface, surface
plt.ylabel('Radial Distibution Function')
plt.xlabel('Radius from Center (A)')
plt.savefig('rdf_plot.pdf')
plt.show()