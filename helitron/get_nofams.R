library(plyr)

w=read.table('W22__Ver12.HelitronScanner.8080.out', header=F, sep='\t')

allw=read.table('W22__Ver12.HelitronScanner.tabnames.terminal30bp.fa.fai', header=F)
nofam=allw$V1[!allw$V1 %in% w$V1]
write(as.character(nofam), 'W22.Helitron.noB73v4Fam.txt', ncolumns=1)
