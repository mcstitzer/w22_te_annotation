library(rtracklayer)
library(stringr)

## switch to read in from command line?
GENOMENAME='W22__Ver12'
GFFNAME=paste(GENOMENAME, '.out.gff', sep='')

### can be generalized to DNA TEs
TSDLEN=5  ## 5 for LTR TEs

a=import.gff3(GFFNAME) ## switch to read in from command line
### get rid of matches < TSD from start of sequence.
### bedtools will do just fine with those that go past the end of the sequence, but will fail with those before the start
a=a[start(a)>TSDLEN,]
a=sort(a, ignore.strand=T)
up=flank(a, width=TSDLEN, start=T)
down=flank(a, width=TSDLEN, start=F)


export.gff3(up, paste(GENOMENAME, '.up.gff3', sep=''))
export.gff3(down, paste(GENOMENAME, '.down.gff3', sep=''))


