library(rtracklayer)
library(data.table)
library(stringr)
library(dplyr)

ltr=import.gff3('ltr/W22__Ver12.allFamilies.100kbFilter.gff3')
ltr$Name=NULL
tir=import.gff3('tir/w22_tir_copies.withdetectMITE.gff3')
tir$ID=tir$Name
tir$Name=NULL
sine=import.gff3('sine/W22.RST.gff3')
hel=import.gff3('helitron/W22.DHH.gff3')
line=import.gff3('line/W22.LINE.gff3')  ### this one is not finished!!

te=c(ltr, tir, sine, hel, line)
## remove identical hits
## but there aren't any now that I've filtered things appropriately...
#te=te[-queryHits(findOverlaps(te, drop.self=T, ignore.strand=T, type='equal')[c(F,T)]),]
queryHits(findOverlaps(te, drop.self=T, ignore.strand=T, type='equal'))
          
          
te=sort(sortSeqlevels(te))
te$sup=substr(te$ID,1,3)


## here, I want TEs where either the helitron is entirely within a diff TE, or the other TE is entirely within a helitron
##   so need to switch column names in the as.data.frame call
helNO=anti_join(rbind(data.frame(findOverlaps(te[te$sup != 'DHH',], te[te$sup=='DHH',], ignore.strand=T, type='start')), 
                data.frame(findOverlaps(te[te$sup != 'DHH',], te[te$sup=='DHH',], ignore.strand=T, type='end')), 
                data.frame(findOverlaps(te[te$sup != 'DHH',], te[te$sup=='DHH',], ignore.strand=T, type='any'))), 
              rbind(data.frame(findOverlaps(te[te$sup != 'DHH',], te[te$sup=='DHH',], ignore.strand=T, type='within')),
                setNames(data.frame(findOverlaps(te[te$sup=='DHH',], te[te$sup != 'DHH',], ignore.strand=T, type='within')), c('subjectHits', 'queryHits'))))$subjectHits
slineNO=anti_join(rbind(data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX', 'DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], te[te$sup %in% c('RST','RIT', 'RIL'),], ignore.strand=T, type='start')), 
                  data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX', 'DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], te[te$sup %in% c('RST','RIT', 'RIL'),], ignore.strand=T, type='end')),
                  data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX', 'DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], te[te$sup %in% c('RST','RIT', 'RIL'),], ignore.strand=T, type='any'))),
                rbind(data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX', 'DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], te[te$sup %in% c('RST','RIT', 'RIL'),], ignore.strand=T, type='within')),
                  setNames(data.frame(findOverlaps(te[te$sup %in% c('RST','RIT', 'RIL'),], te[te$sup %in% c('RLC', 'RLG', 'RLX', 'DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], ignore.strand=T, type='within')), c('subjectHits', 'queryHits'))))$subjectHits
tirNO=anti_join(rbind(data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX'),], te[te$sup %in% c('DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], ignore.strand=T, type='start')), 
                data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX'),], te[te$sup %in% c('DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], ignore.strand=T, type='end')),
                data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX'),], te[te$sup %in% c('DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], ignore.strand=T, type='any'))),
              rbind(data.frame(findOverlaps(te[te$sup %in% c('RLC', 'RLG', 'RLX'),], te[te$sup %in% c('DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], ignore.strand=T, type='within')),
                setNames(data.frame(findOverlaps(te[te$sup %in% c('DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'),], te[te$sup %in% c('RLC', 'RLG', 'RLX'),], ignore.strand=T, type='within')), c('subjectHits', 'queryHits'))))$subjectHits
                 
                 
filtTE=te
filtTE=filtTE[-which(filtTE$sup=='DHH')[unique(helNO)],]
filtTE=filtTE[-which(filtTE$sup %in% c('RST', 'RIT', 'RIL'))[slineNO],]         
filtTE=filtTE[-which(filtTE$sup %in% c('DTA', 'DTC', 'DTM', 'DTH', 'DTT', 'DTX'))[tirNO],]                

filtTE$sup=NULL

filtTE=sort(sortSeqlevels(filtTE), ignore.strand=T)

## output the TEs.
export.gff3(filtTE, 'W22.allTE.gff3')

### get ready to disjoin
te$sup=substr(te$ID, 1,3)  ## get the superfamily three letter code
te.dj=disjoin(te, ignore.strand=T)    ### disjoin overlapping TEs so each bp in genome is only covered once.
## relate these disjoined bits to the TE they belong to. 
te.o=findOverlaps(te.dj, te, ignore.strand=T)
mcols(te.dj)=splitAsList(mcols(te)$ID[subjectHits(te.o)],queryHits(te.o))
mcols(te.dj)$ID=sapply(mcols(te.dj)$X, function(x) paste(as.character(x), collapse=','))
mcols(te.dj)$ID=sapply(te.dj$ID, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])
te.dj$sup=substr(te.dj$ID, 1,3)
## complete helitrons
th=findOverlaps(te.dj, te[te$sup=='DHH',], type='equal')
te.dj$hel=F
te.dj$hel[queryHits(th)]=T
## complete tirs
tt=findOverlaps(te.dj, te[substr(te$sup, 1,2)=='DT',], type='equal')
te.dj$tir=F
te.dj$tir[queryHits(tt)]=T
## complete line/sine
tls=findOverlaps(te.dj, te[te$sup %in% c('RIT', 'RIL', 'RST'),], type='equal')
te.dj$sline=F
te.dj$sline[queryHits(tls)]=T
## complete ltr
te.dj$ltr=F
te.dj$ltr[sapply(1:length(te.dj), function(x) sum(te.dj$hel[x], te.dj$tir[x], te.dj$sline[x]))==0]=T
## intact ltr
ti=findOverlaps(te.dj, te[te$sup %in% c('RLC', 'RLG', 'RLX'),], type='equal')
te.dj$intactltr=F
te.dj$intactltr[queryHits(ti)]=T
## solo ltr
ts=findOverlaps(te.dj, te[te$type=='solo_LTR',], type='any')
te.dj$solo=F
te.dj$solo[queryHits(ts)]=T
thi=findOverlaps(te.dj, te[te$sup=='DHH',], type='equal')
te.dj$intacthel=F
te.dj$intacthel[queryHits(thi)]=T
tti=findOverlaps(te.dj, te[substr(te$sup, 1,2)=='DT',], type='equal')
te.dj$intacttir=F
te.dj$intacttir[queryHits(tti)]=T

## more useful is any that are intact
tei=findOverlaps(te.dj, te, type='equal')
te.dj$intact=F
te.dj$intact[queryHits(tei)]=T

## confirm no overlaps between helitron and others
findOverlaps(te.dj[te.dj$hel,], te.dj[!te.dj$hel,])


                 
export.gff3(te.dj, 'W22.allTE.disjoined.gff3')                 
                 
                 
