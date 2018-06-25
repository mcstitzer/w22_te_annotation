import sys
import textwrap

def readFasta(filename):
        fastafile=open(filename, 'r')
        fastadict={}
        for line in fastafile:          
                if line.startswith('>') and line.strip()[1:] not in fastadict:
                        seqname=line.strip()[1:]
                        fastadict[seqname]=[]
                elif line.startswith('>') and line.strip()[1:] in fastadict:
                        seqname=line.strip()[1:]
#                       continue
                else:
                        fastadict[seqname].append(line.strip())
        for entry in fastadict:
                fastadict[entry]=''.join(fastadict[entry])
        return fastadict

def printFasta(sequence, width=70):
        return '\n'.join(sequence[i:i+width] for i in range(0,len(sequence), width))

tes=open(sys.argv[2], 'r')
te={}
for line in tes:
	fields=line.strip().split('\t')
	te[fields[1]]=fields[0]
tes.close()

#print te

b=readFasta(sys.argv[1])
#out=open('stupid_test.txt', 'w')
#out.write('read in all the fasta entries')
for entry in b: 
	print '>'+te[entry]+'\n'+printFasta(''.join(b[entry]))



