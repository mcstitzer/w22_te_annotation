library(rtracklayer)
library(plyr)

tf=readRDS('../../tf_granges_intact.RDS')
ltr=tf[tf$source=='LTRharvest',]

dig=import.gff3('B73V4.pseudomolecule.ltrdigest.ALL.gff3')

## figure out which entries are LTR retrotransposon entries (so should have same start and end coordinates as ltrharvest)
dig$LTR=substr(dig$ID, 1,3)=='LTR'
dig$LTR[is.na(dig$LTR)]=F
dig$MCSte[dig$LTR]=mapvalues(paste(seqnames(dig), start(dig), end(dig))[dig$LTR], from=paste(seqnames(ltr), start(ltr), end(ltr)), to=as.character(ltr$ID), warn_missing=F)

## get the MCS names for the LTR copies
dig$mcs=mapvalues(as.character(dig$Parent), from=dig$ID[dig$LTR], to=dig$MCSte[dig$LTR], warn_missing=F)
## switch repeat regions to the right MCS name
dig$mcs[substr(dig$mcs, 1,13)=='repeat_region' & !is.na(dig$mcs)]=mapvalues(dig$mcs[substr(dig$mcs, 1,13)=='repeat_region' & !is.na(dig$mcs)], from=dig$mcs[!is.na(dig$MCSte)], to=dig$MCSte[!is.na(dig$MCSte)], warn_missing=F)
## switch MCS names for the final repeat region entry
dig$mcs[is.na(dig$mcs)]=dig$MCSte[!is.na(dig$MCSte)]

digout=dig
digout$MCSte=NULL
digout$LTR=NULL
export.gff3(digout, 'B73V4.pseudomolecule.ltrdigest.MCSnames.gff3')

alllist=list(digout)

for (i in 1:80){
  d=import.gff3(paste('B73V4.pseudomolecule.subtract', i, '.ltrdigest.contignames.gff3.contigpositions.gff3', sep=''))
  d$LTR=substr(d$ID, 1,3)=='LTR'
  d$LTR[is.na(d$LTR)]=F
  d$MCSte[d$LTR]=mapvalues(paste(seqnames(d), start(d), end(d))[d$LTR], from=paste(seqnames(ltr), start(ltr), end(ltr)), to=as.character(ltr$ID), warn_missing=F)
## get the MCS names for the LTR copies
  d$mcs=mapvalues(as.character(d$Parent), from=d$ID[d$LTR], to=d$MCSte[d$LTR], warn_missing=F)
## switch repeat regions to the right MCS name
  d$mcs[substr(d$mcs, 1,13)=='repeat_region' & !is.na(d$mcs)]=mapvalues(d$mcs[substr(d$mcs, 1,13)=='repeat_region' & !is.na(d$mcs)], from=d$mcs[!is.na(d$MCSte)], to=d$MCSte[!is.na(d$MCSte)], warn_missing=F)
## switch MCS names for the final repeat region entry
  d$mcs[is.na(d$mcs)]=d$MCSte[!is.na(d$MCSte)]
  dout=dig
  dout$MCSte=NULL
  dout$LTR=NULL
#  export.gff3(digout, 'B73V4.pseudomolecule.ltrdigest.MCSnames.gff3')
  alllist[[i+1]] = dout
}  

all=do.call(c, alllist)
export.gff3(all, 'B73V4.pseudomolecule.ltrdigest.ALLsubtracted.MCSnames.gff3')
