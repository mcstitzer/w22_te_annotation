import sys

fastatoswitch=open(sys.argv[1], 'r')

newnames=open(sys.argv[2], 'r')


## make a dictionary with old name as key, new name as value
switchdict={}
for line in newnames:
	fields=line.strip().split('\t')
	switchdict[fields[0]]=fields[1]

for line in fastatoswitch:
	if line.startswith('>'):
		print ">"+switchdict[line.strip()[1:]]
	else:
		print line.strip()


