import numpy as np
import sys
import matplotlib.pyplot as plt
import numpy.linalg as npl
from qd_helper import *
import copy
from geom_helper import *
from pdos_helper import dos_grid_general,get_alpha,get_ao_ind, get_ind_ao_underc

#the goal of this code is to make a non-janky way of assessing the different atomic contributions to the DOS of each atom in a given MO
#currently only handles In, Ga, P, F, Cl
#nore that dos contributions from ligands are ignored

QD_file=sys.argv[1] #your .xyz file
bas_file = sys.argv[2]      # txt file with total number of orbitals and occupied orbitals
coeff_file = sys.argv[3]    # txt version of qchem 53.0 OR numpy version
target = sys.argv[4] #integer index of MO you want to decompose

coords,atoms = read_input_xyz(QD_file)

nbas,nocc=np.loadtxt(bas_file,dtype=int,unpack=True)
homo=nocc-1

orb_per_atom={'In': 26, 'Ga': 32, 'P': 18, 'Cl': 18, 'Br': 32,'F':14,'H':5,'O':14, 'C':14,'S':18, 'Li':9,"Al":18,"Zn":31,"S":18,"Se":32,"Si":18}

coeff_file_expand=np.load(coeff_file)
mo_mat=coeff_file_expand['arr_0'] # normalized
mo_e = coeff_file_expand['arr_1']

ind_In = atoms=='In'
ind_P = np.logical_or((atoms=='P'), np.logical_or((atoms=='S'),(atoms=='Se')))
ind_Ga = atoms=="Ga"
ind_Al = np.logical_or(atoms=="Zn", atoms=="Al")
ind_cat=np.logical_or(ind_Al,np.logical_or(ind_In,ind_Ga))
ind_lig = np.logical_or(atoms == "F",atoms=="Cl")
ind_lig2 = np.logical_or(np.logical_or(atoms=="Si",atoms=="O"),atoms=="C")

ga = False
if np.any(ind_Ga):
	ga = True

al = False
if np.any(ind_Al):
	al = True

#math
ind_in_ao,ind_p_ao,ind_ga_ao,ind_al_ao,ind_lig_ao,ind_lig2_ao = get_ao_ind([ind_In,ind_P,ind_Ga,ind_Al,ind_lig,ind_lig2],atoms,nbas,orb_per_atom)
'''
get squared coefficients on core,shell
'''
alpha_list = get_alpha(mo_mat,[ind_in_ao,ind_p_ao,ind_ga_ao,ind_al_ao,ind_lig_ao,ind_lig2_ao])
alpha_in,alpha_p,alpha_ga,alpha_al,alpha_lig,alpha_lig2=alpha_list
alpha_inga=alpha_in+alpha_ga

# make sure alphas add to 1
test = np.all(np.isclose(alpha_in+alpha_p+alpha_ga+alpha_al+alpha_lig+alpha_lig2,1))
if test == False: raise ValueError('Alpha doesnt add to 1!')
print('Alphas add to 1?:',test)
'''
calculate projected DOS
'''
mo_e = mo_e * 27.2114 # MO energy, in eV
E_grid = np.arange(-50,50,0.01) # energy grid to evaluate the DOS over
sigma=0.1 # broadening parameter
print('HOMO energy:',mo_e[homo])
print('HOMO-LUMO gap:', round(mo_e[homo+1]-mo_e[homo],3), "eV")
dos_list = dos_grid_general(E_grid, sigma,mo_e, alpha_list)
in_dos,p_dos,ga_dos,al_dos,lig_dos,lig2_dos=dos_list
inga_dos=in_dos+ga_dos


#some specific code for one dot 
#core_p4c_total=0
#facial_p4c_total=0
#edge_p4c_total=0
#p_3c_total=0
#in3c_oc=0
#in3c_oe=0
#in4c_ic=0
#in_4c_ie=0
#in_core=0
#p_special=0

in_count=0
p_count=0
ga_count=0
al_count=0
s_count=0
for i,atom in enumerate(atoms):
	if ind_cat[i]:
	#if ind_P[i]:
	#if not ind_lig[i]:
	#if ind_lig2[i]:
		in_index_uc=np.zeros(len(atoms),dtype=bool)
		p_index_uc=np.zeros(len(atoms),dtype=bool)

		if ga:
			ga_index_uc=np.zeros(len(atoms),dtype=bool)
		if al:
			al_index_uc=np.zeros(len(atoms),dtype=bool)


		if atom == "In":
			in_index_uc[i]=True
			in_count=in_count+1
		if ind_P[i]:
			p_index_uc[i]=True
			if atom=="P":
				p_count=p_count+1
			if atom=="S" or atom=="Se":
				s_count=s_count+1
		if ga:
			if ind_Ga[i]:
				ga_index_uc[i]=True
				ga_count=ga_count+1
		if al:
			if ind_Al[i]:
				al_index_uc[i]=True
				al_count=al_count+1


		if al:
			al_uc_ind_ao, p_uc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,al_index_uc,p_index_uc)
		if ga:
			ga_uc_ind_ao, p_uc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,ga_index_uc,p_index_uc)
		in_uc_ind_ao,p_uc_ind_ao = get_ind_ao_underc(atoms,nbas,orb_per_atom,in_index_uc,p_index_uc)			

		p_underc = p_uc_ind_ao #np.logical_or(sshell_underc_ind_ao,sshell_underc_ind_amb_ao)
		in_underc = in_uc_ind_ao #np.logical_or(cdshell_underc_ind_ao,cdshell_underc_ind_amb_ao)
		if ga:
			ga_underc = ga_uc_ind_ao
		if al:
			al_underc = al_uc_ind_ao

		p_nouc = np.logical_xor(p_underc,ind_p_ao)
		in_nouc = np.logical_xor(in_underc,ind_in_ao)
		if ga:
			ga_nouc = np.logical_xor(ga_underc,ind_ga_ao)
		if al:
			al_nouc = np.logical_xor(al_underc,ind_al_ao)

		if al:
			alpha_in_uc,alpha_p_uc,alpha_ga_uc,alpha_al_uc,alpha_in_fc,alpha_p_fc,alpha_ga_fc,alpha_al_fc = get_alpha(mo_mat,[in_underc,p_underc,ga_underc,al_underc,in_nouc,p_nouc,ga_nouc,al_nouc])
			alpha_fc=alpha_in_fc+alpha_p_fc+alpha_ga_fc+alpha_al_fc
		elif ga:
			alpha_in_uc,alpha_p_uc,alpha_ga_uc,alpha_in_fc,alpha_p_fc,alpha_ga_fc = get_alpha(mo_mat,[in_underc,p_underc,ga_underc,in_nouc,p_nouc,ga_nouc])
			alpha_fc=alpha_in_fc+alpha_p_fc+alpha_ga_fc
		else:
			alpha_in_uc,alpha_p_uc,alpha_in_fc,alpha_p_fc = get_alpha(mo_mat,[in_underc,p_underc,in_nouc,p_nouc])
			alpha_fc=alpha_in_fc+alpha_p_fc
		

		#GET WHAT WE WANT
		if atom == "In":
			if float(alpha_in_uc[int(target)-1])>0.02:
				print("In MO",target,",",atom,in_count,"(total index",str(i)+")","has a DOS contribution of",round(float(alpha_in_uc[int(target)-1]),3))
		
			#if in_count==5 or in_count==7 or in_count==11 or in_count==1 or in_count==2 or in_count==4 or in_count==13 or in_count==14 or in_count==29 or in_count==26 or in_count==28 or in_count==31:
			#	in3c_oc=in3c_oc+float(alpha_in_uc[int(target)-1])
			#if  in_count==18 or in_count==27 or in_count==25 or in_count==15 or in_count==30:
			#	in3c_oe=in3c_oe+float(alpha_in_uc[int(target)-1])
			#if in_count==3 or in_count==6 or in_count==12 or in_count==9 or in_count==10 or in_count==22 or in_count==19 or in_count==20 or in_count==24 or in_count==16 or in_count==17 or in_count==23:
			#	in4c_ic=in4c_ic+float(alpha_in_uc[int(target)-1])
			#if in_count==21:
			#	in_core=in_core+float(alpha_in_uc[int(target)-1])



		if ind_P[i]:
			if float(alpha_p_uc[int(target)-1])>0.02:
				if atom=="P":
					print("In MO",target,",",atom,p_count,"(total index",str(i)+")","has a DOS contribution of",round(float(alpha_p_uc[int(target)-1]),3))
				if atom=="S" or atom=="Se":
					print("In MO",target,",",atom,s_count,"(total index",str(i)+")","has a DOS contribution of",round(float(alpha_p_uc[int(target)-1]),3))
			#again, very specific code for specific analysis on specific dots
			#if   p_count==7 or p_count==15 or p_count==17 or p_count==8 or p_count==23 or p_count==27 or p_count==10 or p_count==34 or p_count==35 or p_count==5 or p_count==12 or p_count==18 or p_count==4 or p_count==20 or p_count==24 or p_count==2 or p_count==28 or p_count==29:
			#	edge_p4c_total=edge_p4c_total+float(alpha_p_uc[int(target)-1])
			#if   p_count==13 or p_count==22 or p_count==31 or p_count==19 or p_count==26 or p_count==33 or p_count==14 or p_count==21 or p_count==32 or p_count==16 or p_count==25 or p_count==30: 
			#	facial_p4c_total=facial_p4c_total+float(alpha_p_uc[int(target)-1])
			#if p_count==3 or p_count==9 or p_count==11 or p_count==1:
			#	p_3c_total=p_3c_total+float(alpha_p_uc[int(target)-1])


		if ga:
			if ind_Ga[i]:
				if float(alpha_ga_uc[int(target)-1])>0.02:
					print("In MO",target,",",atom,ga_count,"(total index",str(i)+")","has a DOS contribution of",round(float(alpha_ga_uc[int(target)-1]),3))

			#again, very specific code for specific analysis on specific dots
			#if  ga_count==9 or ga_count==51 or ga_count==27 or ga_count==40 or ga_count==12 or ga_count==39 or ga_count==23 or ga_count==26 or ga_count==14 or ga_count==2 or ga_count==44 or ga_count==31:
			#	in3c_oc=in3c_oc+float(alpha_ga_uc[int(target)-1])
			#if  ga_count==4 or ga_count==47 or ga_count==52 or ga_count==38 or ga_count==32 or ga_count==3 or ga_count==1 or ga_count==15 or ga_count==21 or ga_count==50 or ga_count==20 or ga_count==35:
			#	in3c_oe=in3c_oe+float(alpha_ga_uc[int(target)-1])
			#if ga_count==24 or ga_count==22 or ga_count==16 or ga_count==41 or ga_count==46 or ga_count==43 or ga_count==28 or ga_count==30 or ga_count==36 or ga_count==8 or ga_count==13 or ga_count==11:
			#	in4c_ic=in4c_ic+float(alpha_ga_uc[int(target)-1])
			#if ga_count==37 or ga_count==25 or ga_count==49 or ga_count==10 or ga_count==45 or ga_count==33 or ga_count==18 or ga_count==7 or ga_count==42 or ga_count==17 or ga_count==6 or ga_count==46:
			#	in_4c_ie=in_4c_ie+float(alpha_ga_uc[int(target)-1])
			#if ga_count==5 or ga_count==34 or ga_count==19 or ga_count==48:
			#	in_core=in_core+float(alpha_ga_uc[int(target)-1])

		if al:
			if ind_Al[i]:
				if float(alpha_al_uc[int(target)-1])>0.02:
					print("In MO",target,",",atom,al_count,"(total index",str(i)+")","has a DOS contribution of",round(float(alpha_al_uc[int(target)-1]),3))


print("total ligand dos is", alpha_lig[int(target)-1])

#for that extra, dot-specific analysis code
print("")
print(round(mo_e[int(target)-1],3))
print("")
#print("Edge P-4c",round(edge_p4c_total*100,1), "%,",round(100*edge_p4c_total/18,1), "% av")
#print("Face P-4c",round(facial_p4c_total*100,1), "%,",round(100*facial_p4c_total/12,1), "% av")
#print("P-3c",round(p_3c_total*100,1), "%,",round(100*p_3c_total/4,1), "% av")
#print("P total is",round(p_3c_total*100,3)+round(facial_p4c_total*100,3)+round(edge_p4c_total*100,3)+round(p_special*100,3), "percent")
#print("")
#print("OC Ga-3c",round(in3c_oc*100,1), "%,",round(100*in3c_oc/12,1), "% av")
#print("OE Ga-3c",round(in3c_oe*100,1), "%,",round(100*in3c_oe/12,1), "% av")
#print("IC Ga-4c",round(in4c_ic*100,1), "%,",round(100*in4c_ic/12,1), "% av")
#print("IE Ga-4c",round(in_4c_ie*100,1), "%,",round(100*in_4c_ie/12,1), "% av")
#print("Core Ga-4c",round(100*in_core,1), "%,",round(100*in_core/4,1), "% av")