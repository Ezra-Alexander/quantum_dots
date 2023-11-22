#!/bin/bash

#The laptop portion of the script that interpolates cut-outs of any distorted atom from one dot to another
#can handle a variable number of atoms
#the script takes the target atom given to it plus everything bonded to those atoms in both dots (ideally the same number)
#meant to work where the "pristine" is the bulk crystal (P can be substituted to ligands as needed) and the "defect" is the final structure with the distorted atoms
#can also work where the initial structure is another dot (say, before defect)
#in fact, the bulk crystal initial state only really works if i'm doing a single target atom
#currently only works for In, Ga, P. Anything else in the structure is indexed as a ligand. It wouldn't be hard to fix tho
#this outputs a single .xyz. On the cluster, use interpolation_runner.sh to slice it up and run the single points

#update:
#P-terminated cutouts will by default be negatively charged. Everything should be Li terminated

#update:
#geodesic-interpolator fails if the atoms in the .xyz files are not in the same order in the initial and final state

#update:
#the atoms in the initial frame need to be in the same order as the closest corresponding atom in the final frame

pristine=$1 #.xy
defect=$2 #.xyz
target=$3 #the species you are targeting
#cation=$3
#then supply the target atoms as arguments (atom-specific indeces), accepts multiple

shift 3 #all subsequent arguments are the atom-specific indeces of atoms to cut out



#python3 ~/wormk/code/p_carver.py $pristine $defect "$@" #outputs endpoints.xyz, but these still are P-In(Ga), not P-Li

# python3 ~/wormk/code/general_carver.py $pristine $defect $target "$@"

python3 ~/wormk/code/in3c_carver.py $pristine $target "$@" #for this one, the first index is the initial atom and the second is the final atom, both in the pristine .xyz (defect .xyz if fully skipped)


#comment this out for In/Ga targets
# if [ "$target" = "P" ]
# then
# 	python3 ~/wormk/code/bond_adjuster_p4c.py endpoints.xyz endpoints.xyz In P 2.265 #makes them P-Li. 2.384 is the P-Li bond length I was using, but i calculated 2.265 to be more optimal
# fi

geodesic_interpolate endpoints.xyz