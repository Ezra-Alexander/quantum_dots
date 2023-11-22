#!/bin/bash

# a nice simple script that takes a single point input and runs Pipek-Mezey localizations on specified orbitals
# it would be nice if the orbitals could be specified in the "LUMO+X" convention instead of absolute number
# the way this code is currently written, it will not mix occupieds with virtuals (which is probably fine)

#again, this runs only a single job, and should be called from within the appropriate "localization" directory
#is currently configured to always run pipek-mezey localizations with the rij-def2-svp auxilliary basis set
#now, you might think I don't need the submit script because I don't need another PDOS, but its nice to have for dos component analysis

#the qchem .in for the single point plot you already ran
input=$1

#good old nbas.txt
nbas=$2

#the name of the single point scratch directory. should be 2 directories back. Don't copy it, and don't include the dots
scr=$3

#whether the states of interest are occupied (o) or unoccupied (u)
occ=$4

#the relative index of all states that you want to localize are then specified. Can be any number of states. 0 for HOMO/LUMO
shift 4

#make the localize.in and the submit script
python3 ~/code/localizer.py $input $nbas $scr $occ "$@"

#run the job
sbatch submit_loc.sh