import sys
import numpy
import math

file = sys.argv[1]

to_average = []
#average = 0
#count = 0
with open(file,'r') as file:
	for i,line in enumerate(file):
		# for DMEs
		if i > 1:
			if int(line.strip().split()[0]) != int(line.strip().split()[3]):
				#configured for MT pristine
				if (int(line.strip().split()[0])<1752 and int(line.strip().split()[3])<1752) or (int(line.strip().split()[0])>1751 and int(line.strip().split()[3])>1751):
					to_average.append(float(line.strip().split()[12]))
				#configured for MC pristine
				#if (int(line.strip().split()[0])<1848 and int(line.strip().split()[3])<1848) or (int(line.strip().split()[0])>1847 and int(line.strip().split()[3])>1837):
				#	to_average.append(float(line.strip().split()[12]))
		#for IPRs
		#if float(line.strip().split()[0]) < 1:
		#	to_average.append(float(line.strip().split()[0]))
		#average = average + float(line.strip().split()[0])
		#count = count+1
	#average = average / count

print(numpy.mean(to_average))
#print(average)
#print(max(to_average))