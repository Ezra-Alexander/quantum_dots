import numpy as np
import sys

#print HOMO energy and HOMO-LUMO Gao, and the dipole matrix elements (squared) between the trap and pre-selected band edge states
#as written this only works for one trap, but can do LUMO+1/HOMO-1
#low_orb_ipr.csv should be in the directory

#the aims .out
out = sys.argv[1]

#the aims dipole matrix
dipole = sys.argv[2]

#The lowest OMO to get: 0 for HOMO, 1 for HOMO-1, etc
omo = sys.argv[3]

#the highest UMO to get: 0 for LUMO, 1 for LUMO+1, etc
umo = sys.argv[4]

#what kind of trap we have. If the trap is the homo should be "o"; if the trap is the lumo should be "u"; if the trap is lumo+1 it should be "u2"; homo-1 "o2"
trap = sys.argv[5]

#trap index in def2-svp (different than trap index in Aims)
q_ti = int(sys.argv[6])


#step one - go through .out
#want to print HOMO energy and HOMO-LUMO gap
#also want to save the HOMO index for future use

hl_gap = 0
homo_e = 0
t_index = 0
flag = 0
eig_flag = False
eig = []
converged_flag = False
with open(out, "r") as o:
	for i, line in enumerate(o):
		if line.find('Overall HOMO-LUMO gap:') != -1:
			hl_gap = float(line.strip().split()[3])
		if line.find('Highest occupied state (VBM) at') != -1:
			homo_e = float(line.strip().split()[5])
			flag = 0
		if line.find("Writing Kohn-Sham eigenvalues.") != -1 and flag == 0:
			flag = 1
		if line.find("Self-consistency cycle converged.") != -1:
			converged_flag = True
		if eig_flag and line =="\n":
			eig_flag = False
		if eig_flag:
			eig.append(float(line.strip().split()[3]))
		if converged_flag and line.find("State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]") != -1:
			eig_flag = True
		if line.find("2.00000") != -1 and flag == 1:
			flag = 2
		elif flag == 2:
			if line.find("0.00000") != -1:
				t_index = int(line.strip().split()[0])
				flag = 1
			elif line.find("2.00000") != -1:
				flag = 2
			else: 
				flag = 1

if trap == "o":
	t_index = t_index-1
elif trap =="u2":
	t_index = t_index+1
elif trap=="o2":
	t_index=t_index-2
elif trap=="u3":
	t_index=t_index+2

print()
print("HOMO energy is ", homo_e, " eV")
print()
print("HOMO-LUMO gap is ", hl_gap, " eV")
print()

lumo_e = int(homo_e)+int(hl_gap)


#step 2 - extract all the dipole matrix elements we want


t_els = []
with open(dipole, "r") as d:
	for i, line in enumerate(d):
		if line.strip().split()[0] == str(t_index) or line.strip().split()[3] == str(t_index):
			t_els.append(line.strip().split())

just_d = []
for i, line in enumerate(t_els):
	just_d.append(line)


lb=0
ub=0
if trap == "o":

	if int(omo)==0:
		print("Wrong OMO bound!")

	lb = t_index - int(omo)
	ub = t_index + int(umo) + 1

if trap == "u":

	if int(umo)==0:
		print("Wrong UMO bound!")

	lb = t_index - int(omo) - 1
	ub = t_index + int(umo)

if trap == "u2":
	if int(umo)==0 or int(umo)==1:
		print("Wrong UMO bound!")

	lb = t_index - int(omo) - 2
	ub = t_index + int(umo) - 1

if trap == "u3":
	if int(umo)==0 or int(umo)==1 or int(umo)==2:
		print("Wrong UMO bound!")

	lb = t_index - int(omo) - 3
	ub = t_index + int(umo) - 2


if trap == "o2":
	if int(omo)==0 or int(omo)==1:
		print("Wrong OMO bound!")

	lb = t_index - int(omo) + 1
	ub = t_index + int(umo) + 2


test = False
while test == False:

	count = 0

	for i,line in enumerate(t_els):
		if int(line[0]) < lb or int(line[3]) > ub or int(line[0]) == int(line[3]):
			t_els.pop(i)
		else:
			count = count + 1

	if count == len(t_els):
		test = True

print("Our trap is MO number", t_index)
print()

for i, line in enumerate(t_els):
	print("The dipole matrix element ((x_ij^2+y_ij^2+z_ij^2)/3) between", line[0], "and", line[3], "is", line[12])
	print()


#step 3: open low_orb_ipr and get the relevant IPRs, determine appropriate thresholds and which MOs pass

with open("low_orb_ipr.csv", "r") as ipr:
	lines = ipr.readlines()
	band_edge = []

	if trap == "u":
		band_edge = lines[q_ti-1:q_ti+50]

	if trap == "u2":
		band_edge = lines[q_ti-2:q_ti+49]

	if trap == "u3":
		band_edge = lines[q_ti-3:q_ti+48]

	if trap == "o":
		band_edge = lines[q_ti-51:q_ti]
		band_edge.reverse()

	if trap == "o2":
		band_edge = lines[q_ti-50:q_ti+1]
		band_edge.reverse()

	for i, line in enumerate(band_edge):
		band_edge[i] = float(band_edge[i])

#	cutoff = 0
#	for i, line in enumerate(band_edge):
#		if i > 1:
#			if line/band_edge[i-1] >= 2:
#				cutoff = i
#				break

#	threshold = (np.mean(band_edge[:cutoff])+1.75*np.mean(band_edge[cutoff:cutoff+3]))/2

#let's just use general thresholds for the conduction and valence bands based on atom populations

with open("geometry.in", "r") as geom:
	cation_count = 0
	anion_count = 0
	total_count = 0
	for i,line in enumerate(geom):
		if line.find(" P") != -1 and line.find("atom") != -1:
			anion_count = anion_count+1
			total_count=total_count+1
		elif line.find(" In") != -1 and line.find("atom") != -1:
			cation_count=cation_count+1
			total_count=total_count+1
		elif line.find("atom") != -1:
			total_count = total_count+1

	thresh = 0
	if trap == "o" or trap =="o2":
		thresh = 0.5*anion_count/total_count

	if trap == "u" or trap == "u2" or trap == "u3":
		thresh = 0.4*cation_count/total_count

	#this all simply doesn't work for the magic cluster, it has an average IPR ~= twice MS
	#gotta comment this out for other structures
	#thresh = 0.16


	print("IPR threshold is:", thresh)
	print()

#it would be nice to have 1 array that has all the info for each MO

mos = []
count = t_index

if trap == "o2":
	count = count + 1

if trap == "u2":
	count = count - 1

if trap == "u3":
	count = count - 2


for i, line in enumerate(band_edge):

	d_el = 0
	for j, el in enumerate(just_d):
		if count == t_index:
			if el[0] == el[3]:
				d_el = el[12]
		else:
			if int(el[0]) == count or int(el[3]) == count:
				d_el = el[12]



	mos.append([count, line, eig[count-1], float(d_el)])

	if trap == "o2" or trap == "o":
		count = count - 1

	if trap == "u2" or trap == "u" or trap == "u3":
		count = count + 1

print(mos[0])
print()

#step 4: calculate and print key metrics (low_energy and high_energy averages above certain IPR thresholds)

t_e = 0
for i, mo in enumerate(mos):
	if int(mo[0]) == int(t_index):
		t_e = float(mo[2])

#test first energy metric
test = False
e_range = 0.5
while test==False:

	count = 0

	tested = []
	for i, mo in enumerate(mos):
		if abs(float(mo[2]) - t_e) < e_range and float(mo[1]) >= thresh:  #and int(mo[0])>t_index: #last condition only when you want to skip LUMO
			count = count +1 #count cound include the trap and MOs on either side of it
			tested.append(float(mo[3]))

	#print(t_e+e_range)

	if count > 2:
		test = True
	else:
		e_range = e_range+0.25



print("The average DME between the trap and all MOs within", e_range, "eV with IPR >", round(thresh,3), "( of which there are", count, ") is ", round(np.mean(tested),3))
print()

count = 0
e_range = e_range+0.5
#print(t_e+e_range)
test2 = []
for i, mo in enumerate(mos):
		if abs(float(mo[2]) - t_e) < e_range and float(mo[1]) >= thresh: #and int(mo[0])>t_index: #last condition only when you want to skip LUMO
			count = count +1 #count cound include the trap and MOs on either side of it
			print(mo)
			print()
			test2.append(float(mo[3]))

print("The average DME between the trap and all MOs within", e_range, "eV with IPR >", round(thresh,3), "( of which there are", count, ") is ", round(np.mean(test2),3))
print()

