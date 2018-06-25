library(rtracklayer)
library(stringr)
library(data.table)
library(plyr)

GENOMENAME='W22'
#GENOMENAME=commandArgs(trailingOnly=TRUE)
SHORTID='Zm00004b'



######### Read in silix results (family assignments)
## these match to B73 families
a=read.table('W22-matches.noTSD..MCSnames.id80.cov80.out', header=F)

## these are new families in W22
n=read.table('W22.RST.noB73v4fam.8080.fnodes', header=F)

## find the maximum number in B73, increment.
allSine=fread('grep ">" RST.B73v4.fa', header=F)
allSine$famNum=substr(allSine$V1,5,9)
maxFamNum=as.numeric(max(allSine$famNum))

### Switch SINE to RST, the wicker 3 letter code
## because these are plants, and found with the rna pol A and B boxes, these are all RST
## remove one digit because there are not 100,000 families, one order of mag less
#f$V1=sub('SINE0', 'RST', f$V1)  ### let's switch and do as early as we can!!!
#a$V1=sub('SINE0', 'RST', a$V1)


######### Read in sine_finder results, fix naming, classes
sine=read.table(paste(GENOMENAME, '-matches.csv', sep=''), header=T)
sine=unique(sine)   ## concerning - does it append to file? 
sine=sine[-which(sine$name=='name'),]    ## because this looks like a new header line!!
sine$sequence=NULL
sine$direct='+'
sine$namecomp=paste(sine$name, sine$start, sine$end, 'F', paste('TSDlen', sine$TSD.len, sep=''), paste('TSDscore', sine$TSD.score, sep=''), paste('TSDmism', sine$TSD.mism, sep=''), sep='_')

sine$B73fam=mapvalues(sine$namecomp, from=a$V1, to=substr(a$V2,1,8))
sine$B73fam[grepl('TSD', sine$B73fam)]=NA

sine$W22=as.character(mapvalues(sine$namecomp, from=n$V2, to=as.character(n$V1)))
sine$W22[grepl('TSD', sine$W22)]=NA


## first for new w22 families  ## this is super slow!!
for (x in 1:length(table(sine$W22))){
  famNum=str_pad(maxFamNum + x , 5, pad='0') # we want to increment families
  famName=names(rev(sort(table(sine$W22))))[x]  ## and keep track of original families
  sine$w22fam[sine$W22==famName & !is.na(sine$W22)]=paste('RST', famNum, SHORTID, str_pad(1:sum(sine$W22[!is.na(sine$W22)]==famName), 5, pad='0'), sep='')
}              
           
## assign 11 digit copy name (SHORTIDXXXXX) for each copy in an existing B73 family
sine$Name=NA
for (x in names(table(sine$B73fam))){
  sine$Name[sine$B73fam==x & !is.na(sine$B73fam)]=paste(x, SHORTID, str_pad(1:sum(sine$B73fam[!is.na(sine$B73fam)]==x), 5, pad='0'), sep='')
}
sine$Name[is.na(sine$B73fam)]=sine$w22fam[is.na(sine$B73fam)]                   

## done!!!!
sum(is.na(sine$Name))
                                




## output the gff

sine.gff=data.frame(sine$name, 'SineFinder', 'SINE_element', sine$start, sine$end, '.', '+', '.', paste('ID=', sine$Name, sep=''))
write.table(sine.gff, paste(GENOMENAME, '.RST.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

### also keep track of fasta names and newly assigned gffnames as to easily convert between the two (e.g. switching fasta names)
write.table(data.frame(sine$Name, sine$namecomp), paste(GENOMENAME, '.RST.gffname.fastaname.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

