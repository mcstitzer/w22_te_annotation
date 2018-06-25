library(plyr)
source('../CONFIG.R')

w=read.table(paste0(GENOMENAME, '.5ltr.fa.MCSnames.id80.cov80.out'), header=F, sep='\t')


allw=read.table(paste0(GENOMENAME, '.5ltr.fa.fai'), header=F)
nofam=allw$V1[!allw$V1 %in% w$V1]
write(as.character(nofam), paste0(GENOMENAME, '.noB73v4Fam.txt'), ncolumns=1)


#l=read.table('W22.ltrretrotransposon.gff3', header=F)
#l$name=paste(l$V1, l$V4, l$V5, sep='_')

#sum(l$name %in% c(allw$V1, w$V1))
#sum(!l$name %in% c(allw$V1, w$V1))

