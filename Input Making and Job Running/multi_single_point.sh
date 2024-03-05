#!/bin/bash

# a script built on xyz_to_single_point.sh
# but that will loop it over multiple (specified) directories
# and will make the .xyz in the right way from a .out file

#does assume that, for qchem, that the only .out in the directory is the opt.out


cp2k=$1 #'cp2k' if opts were run in cp2k, 'qchem' if opts were run in qchem
descriptor=$2 #a general descriptor string to be added to all single point names
#the following parameters will be applied to all SPs
cores=$3 #the number of cores you want the calculation to run on
plot_min=$4 #the lowest energy MO you want to plot. 0 for HOMO, 1 for HOMO-1, etc
plot_max=$5 #the highest energy MO you want to plot. 0 for LUMO, 1 for LUMO+1, etc
mode=$6 #just a plot ("plot"), just an ipr ("ipr"), or the ipr+plot ("both")
mem_static=$7 #the mem_static value in qchem's $rem, in MB
priority=$8 #the priority for the job. If on Ulysses, say ulysses. priotity automatically sets mem_total

shift 8 #the remaining inputs are the directories to loop over. No /

for var in "$@"
do
    cd "$var"

    name="optzd_${var}_${descriptor}.xyz"

    if [ "$cp2k" = "cp2k" ]
    then    	
    	python3 ~/code/make_cp2k_xyz.py *-pos-1.xyz $name
    elif [ "$cp2k" = "qchem" ]
    then
    	python3 ~/code/output_xyz.py *.out $name
    else
    	echo "1st input should be either 'qchem' or 'cp2k'"
    fi

    sh ~/code/xyz_to_single_point.sh $name $cores $plot_min $plot_max $mode $mem_static $priority

    cd ../
done