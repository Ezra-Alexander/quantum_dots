import sys

#written for ipr.csv but works for anything tbh

file = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

targets = range(start, end+1)
iprs=[]

with open(file,'r') as file:
	for i,line in enumerate(file):
		if i+1 in targets:
			print(str(i+1)+":", round(float(line.strip()),3))
			iprs.append(float(line.strip()))

print()
print("Average IPR over interval:", round(sum(iprs)/len(iprs),3))