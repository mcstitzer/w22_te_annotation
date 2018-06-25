import sys
import re

### takes gff from genometools, writes gff with contig names to stdout

# genome tools renames sequences so it can search by name
##     ##sequence-region   seq0 1 106474
##     #utg10007
##     seq0    LTRharvest      LTR_retrotransposon     58963   74561   .       ?       .       ID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=96.31;seq_number=0

gff=open(sys.argv[1], 'r')

seq=[]
contig=[]


## read in des file

seqindex=0
if len(sys.argv) > 2:
	des=open(sys.argv[2], 'r')
	for line in des:
		seq.append('seq'+str(seqindex))
		contig.append(line.strip())
		seqindex+=1    ## increment up
	des.close()
	linecount=0
#	print contig
#	print seq
	for line in gff:
		fields=line.strip().split()
		if line.startswith('#'):  ## we don't care about these lines because we've got our des file!
			continue
		elif len(fields)==9:
			contigname=contig[seq.index(fields[0])]
			fields[0]=contigname
			print '\t'.join(fields)
		else:
			print line.strip()

elif len(sys.argv)==2:
	linecount=0
	for line in gff:
		fields=line.strip().split()
		if line.startswith('##sequence-region'):
			seq.append(line.split()[1])
		elif line.startswith('#') and len(fields)==3:   ### this is not general!
			contig.append(line.strip().split()[0][1:])
		elif line.startswith('#') and re.search('#\d', line) is not None:
			contig.append(line.strip()[1:])
		elif line.startswith('#chr') or line.startswith('#scaffold'):  ## for ph207, these are the only starts to fasta entries
			contig.append(line.strip()[1:])
		elif line.startswith('#Scrne'): ## for diplo, these are the contigs
			contig.append(line.strip()[1:])
		elif len(fields)==9:
			contigname=contig[seq.index(fields[0])]
			fields[0]=contigname
			print '\t'.join(fields)
		else:
			print line.strip()
gff.close()
