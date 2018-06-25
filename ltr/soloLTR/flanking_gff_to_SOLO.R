library(rtracklayer)
library(stringr)
library(plyr)

GENOME='W22__Ver12'

## generate and up and down sequence
a=import.gff(paste0(GENOME, '.out.gff'))
head(a)
up=flank(a, width=5, start=T)
down=flank(a, width=5, start=F)


## read in the upstream TSD sequence
up.tsd=read.table(paste0(GENOME, '.up.tabout'), header=F)
up.tsd$chr=str_split_fixed(up.tsd$V1, ':', 2)[,1]
up.tsd$start=str_split_fixed(str_split_fixed(up.tsd$V1, ":", 2)[,2], '-', 2)[,1]
up.tsd$end=str_split_fixed(str_split_fixed(up.tsd$V1, ":", 2)[,2], '-', 2)[,2]
up.tsd$start=as.numeric(up.tsd$start)
up.tsd$end=as.numeric(up.tsd$end)
tsdup=GRanges(up.tsd$chr, IRanges(start=up.tsd$start, end=up.tsd$end))
tsdup$tsd=up.tsd$V2
o=findOverlaps(up, tsdup, ignore.strand=T)
up$utsd=NA
up$utsd[queryHits(o)]=as.character(tsdup$tsd[subjectHits(o)])

## and the downstream tsd sequence
down.tsd=read.table(paste0(GENOME, '.down.tabout'), header=F)
down.tsd$chr=str_split_fixed(down.tsd$V1, ':', 2)[,1]
down.tsd$start=str_split_fixed(str_split_fixed(down.tsd$V1, ":", 2)[,2], '-', 2)[,1]
down.tsd$end=str_split_fixed(str_split_fixed(down.tsd$V1, ":", 2)[,2], '-', 2)[,2]
down.tsd$start=as.numeric(down.tsd$start)
down.tsd$end=as.numeric(down.tsd$end)
tsddown=GRanges(down.tsd$chr, IRanges(start=down.tsd$start, end=down.tsd$end))
tsddown$tsd=down.tsd$V2
o=findOverlaps(down, tsddown, ignore.strand=T)
down$dtsd=NA
down$dtsd[queryHits(o)]=as.character(tsddown$tsd[subjectHits(o)])


## merge them with the repeatmasker output
solo=a[up$utsd==down$dtsd & !is.na(up$utsd) & !is.na(down$dtsd),]
solo$tsd=up$utsd[up$utsd==down$dtsd & !is.na(up$utsd) & !is.na(down$dtsd)]
solo$fam=substr(solo$Target, 7,14)
solo$supfam=substr(solo$Target, 7,9)
solo=solo[solo$tsd != 'NNNNN',]

dim(solo)


#### apply 808080
falens=read.table('b73_w22_5ltr.MCSname.fa.fai', stringsAsFactors=F)
falens$V1=gsub("B73v4", "Zm00001d", falens$V1)
solo$Target=gsub("B73v4", "Zm00001d", solo$Target)
solo$ref=substr(solo$Target,7,27)
solo$reflen=as.numeric(mapvalues(solo$ref, from=falens$V1, to=falens$V2, warn_missing=F))

solo$eighty=width(solo)/solo$reflen > 0.8
solo$proportion=width(solo)/solo$reflen
solo=solo[solo$eighty,]

## drop self overlaps, these do not make sense
solo=solo[-queryHits(findOverlaps(solo, drop.self=T, ignore.strand=T)),]

### although some are slightly >1 proportion, we'll consider these to be real, as the solo LTR is slightly bigger than the ref allele.
## write repeatmasker to gff
w=import.gff3('../../w22/W22.allTE.gff3')

## there are some solo LTRs that are the entirety of the TE?? remove them
wi=findOverlaps(w, solo, ignore.strand=T, type='within')
solo=solo[-subjectHits(wi),]

wi=findOverlaps(solo, w, ignore.strand=T, type='within')
solo$wi=countOverlaps(solo, w, ignore.strand=T, type='within')

wo=findOverlaps(solo, w, ignore.strand=T, type='any')

toremove=queryHits(wo)[!subjectHits(wo) %in% subjectHits(wi)]  ## these overlap the edges of another.
solo=solo[-toremove,]

#temps=solo
#mcols(temps)=data.frame(TEID=solo$fam, type='solo')
#tempw=w
#mcols(tempw)=data.frame(TEID=w$ID, type='TE')
#sw=c(temps, tempw)
#swd=disjoin(sw, ignore.strand=T)

## name solos, make format for merge

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
mcols(solo)$Name=NULL
for (x in names(table(mcols(solo)$fam))){
  mcols(solo)$Name[mcols(solo)$fam==x & !is.na(mcols(solo)$fam)]=paste(x, 'Zm00004bS', str_pad(1:sum(mcols(solo)$fam[!is.na(mcols(solo)$fam)]==x), 4, pad='0'),sep='')
}


solos=solo
mcols(solos)=data.frame(source='RepeatMasker', type='solo_LTR', score=NA, phase=NA, ID=solo$Name)

wc=c(w, solos)
wc=sort(wc, ignore.strand=T)


export.gff3(wc, '../../w22/W22.allTE.withSolo.gff3')

