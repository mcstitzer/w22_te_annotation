library(plyr)
w=read.table('W22-matches.noTSD..MCSnames.id80.cov80.out', header=F, sep='\t')


allw=read.table('W22-matches.noTSD.fa.fai', header=F)
nofam=allw$V1[!allw$V1 %in% w$V1]
write(as.character(nofam), 'W22.RST.noB73v4Fam.txt', ncolumns=1)


#l=read.table('W22.ltrretrotransposon.gff3', header=F)
#l$name=paste(l$V1, l$V4, l$V5, sep='_')

#sum(l$name %in% c(allw$V1, w$V1))
#sum(!l$name %in% c(allw$V1, w$V1))

