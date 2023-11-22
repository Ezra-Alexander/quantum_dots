#!/bin/bash

#this bash script should call the titular python script (which edits a qchem optimization input so that it now uses the rca scf algorithm)
#it then submits the edited qchem input to the queue (with 8 cores)
#written for ulysses

input=$1
output=$2
prefix=restart_

python /work/ezraa/code/max_opt_cycle_restartifier.py $input $output

sqthis -J $prefix$input -c 8 -t unlimited qchem.latest -nt 8 $prefix$input $prefix$output
