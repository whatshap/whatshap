import sys

# assumption ... all S' and before L's
filename = sys.argv[1]
#out = sys.argv[2]

d={}
count=1
with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'S'):
			d[var[1]]=count
		count=count+1

with open(filename) as fp:
	for line in fp:
		var=line.split('\t')
		if(var[0] == 'S'):
			print(var[0] + "\t" + str(d[var[1]]) + "\t" + var[2].rstrip())
		else:
			print(var[0] + "\t" + str(d[var[1]]) + "\t" + var[2] + "\t"+ str(d[var[3]]) + '\t'+ var[4] + '\t'+ var[5].rstrip())
