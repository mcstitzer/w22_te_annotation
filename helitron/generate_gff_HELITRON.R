library(rtracklayer)
library(stringr)
library(data.table)
library(plyr)

GENOMENAME='W22'
#GENOMENAME=commandArgs(trailingOnly=TRUE)
SHORTID='Zm00004b'



######### Read in silix results (family assignments)
## these match to B73 families
a=read.table('W22__Ver12.HelitronScanner.8080.out', header=F)

## these are new families in W22
n=read.table('W22.Helitron.noB73v4Fam.8080.fnodes', header=F)

## find the maximum number in B73, increment.
allH=fread('/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/helitron/hel_old_new_out.txt', header=F)
allH$famNum=substr(allH$V2,4,8)
maxFamNum=as.numeric(max(allH$famNum))

### Switch HEL to DHH, the wicker 3 letter code

######### Read in sine_finder results, fix naming, classes
hel=fread('cat W22*.tabout', header=F, sep='\t')

### example line 10      299043  314401  +       10_299043_314401        7       7       TCTATATATACATATTACTCCATGTCTATA  AAGCTTCTCTACCATGGTGGTGCTTGCTAA
names(hel)=c('chr', 'start', 'end', 'orientation', 'ID', 'score5', 'score3', 'fivebp', 'threebp')

hel.gr=GRanges(seqnames=hel$chr, IRanges(start=hel$start, end=hel$end))
## we want to get rid of helitron calls that overlap, because it's difficult to prioritize which (if any) is real.
rmRows=queryHits(findOverlaps(hel.gr, drop.self=T, ignore.strand=T))

hel=hel[-rmRows,]

hel$b73fam=mapvalues(hel$ID, from=a$V1, to=substr(as.character(a$V2),1,8), warn_missing=F)
hel$b73fam[!grepl('DHH', hel$b73fam)]=NA

hel$w22=mapvalues(hel$ID, from=n$V2, to=as.character(n$V1))
hel$w22[!grepl('HEL', hel$w22)]=NA

## first for new w22 families  ## this is super slow!!
for (x in 1:length(table(hel$w22))){
  famNum=str_pad(maxFamNum + x , 5, pad='0') # we want to increment families
  famName=names(rev(sort(table(hel$w22))))[x]  ## and keep track of original families
  hel$w22fam[hel$w22==famName & !is.na(hel$w22)]=paste('DHH', famNum, SHORTID, str_pad(1:sum(hel$w22[!is.na(hel$w22)]==famName), 5, pad='0'), sep='')
}              
           
## assign 11 digit copy name (SHORTIDXXXXX) for each copy in an existing B73 family
hel$Name=NA
for (x in names(table(hel$b73fam))){
  hel$Name[hel$b73fam==x & !is.na(hel$b73fam)]=paste(x, SHORTID, str_pad(1:sum(hel$b73fam[!is.na(hel$b73fam)]==x), 5, pad='0'), sep='')
}
hel$Name[is.na(hel$b73fam)]=hel$w22fam[is.na(hel$b73fam)]                   

## done!!!!
sum(is.na(hel$Name))
                                




## output the gff

hel.gff=data.frame(hel$chr, 'HelitronScanner', 'helitron', hel$start, hel$end, '.', hel$orientation, '.', paste('ID=', hel$Name, sep=''))
write.table(hel.gff, paste(GENOMENAME, '.DHH.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

### also keep track of fasta names and newly assigned gffnames as to easily convert between the two (e.g. switching fasta names)
write.table(data.frame(hel$Name, hel$ID), paste(GENOMENAME, '.DHH.gffname.fastaname.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
