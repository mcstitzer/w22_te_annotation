library(stringr)
library(plyr)
library(rtracklayer)

GENOMENAME='W22__Ver12'
#MAXSUBTRACT=9
MAXSUBTRACT=80
#### read in all the masked TEs

masked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '/', GENOMENAME, '.subtract', 1:MAXSUBTRACT, '.ltrharvest.contignames.gff3.contigpositions.gff3', sep='')
filenames=c(paste('mask_subtract/',GENOMENAME, '.ltrharvest.contignames.gff3', sep=''), masked)

#ltrs=lapply(filenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
ltrs=lapply(filenames, function(x) {a=import.gff(x)  ### no asRangedData in new version of R
									mcols(a)$nest=which(filenames==x)-1
									a=a[mcols(a)$type=='LTR_retrotransposon',]
									mcols(a)$ltr_similarity=as.numeric(mcols(a)$ltr_similarity)
									mcols(a)$nestlevel=which(filenames==x)-1
									return (a)})
ltrs.mask=do.call(c, ltrs)


### add in the coverage of each fragment - how many are overlapping
seqcov=coverage(ltrs.mask)
mcols(ltrs.mask)$cov=NA  ## initialize the variable before loop
for (i in as.character(seqlevels(ltrs.mask))){ ltrs.mask[as.character(seqnames(ltrs.mask))==i,]$cov=as.numeric(seqcov[[i]][start(ltrs.mask[as.character(seqnames(ltrs.mask))==i,])])}

### get the protein results for each
##### MAKE SURE THESE HAVE THE ^B REMOVED!
dignames=paste('mask_subtract/digest_tabout_vim/', GENOMENAME, '.subtract', 1:MAXSUBTRACT, '.ltrdigest_tabout.csv', sep='')
digfiles=c(paste('ltrdigest/', GENOMENAME, '.ltrdigest_tabout.csv', sep=''), dignames)

split.fam=function(s) strsplit(as.character(s), '_')[[1]][1]
prot=lapply(digfiles, function(x) {m=read.table(x, sep='\t', header=T)
									m$sequence=str_split_fixed(as.character(m$sequence), '_', 3)[,1]
									m$rle=''
									print(x)
									m$rle[m$Protein.domain.hits!='']=sapply(which(m$Protein.domain.hits!=''), function(x)  rle(sapply(ldply(strsplit(as.character(m$Protein.domain.hits[x]), '/')), split.fam))$values) 
									m$geneorder=sapply(1:nrow(m), function(x) paste(m$rle[[x]], collapse='.'))
									m$nestlevel=which(digfiles==x)-1
									return(m)
									}
			)
			
prot.mask=do.call(rbind, prot)

## these should be the same order as the ltrs.mask object
## also check that the start values are the same to be super sure
length(ltrs.mask)==nrow(prot.mask)
all(start(ltrs.mask)==prot.mask$element.start) ### this will be false now that these are subtracted!!!!
all(as.character(seqnames(ltrs.mask))==prot.mask$sequence)
mcols(ltrs.mask)$overlaps=countOverlaps(ltrs.mask)
mcols(ltrs.mask)$ownwidth=prot.mask$element.length


ltrs.dj=disjoin(ltrs.mask, ignore.strand=T)
ltrs.overlap=findOverlaps(ltrs.dj, ltrs.mask)
mcols(ltrs.mask)$name=paste(seqnames(ltrs.mask), start(ltrs.mask), end(ltrs.mask), sep='_')
mcols(ltrs.dj)=splitAsList(mcols(ltrs.mask)$name[subjectHits(ltrs.overlap)],queryHits(ltrs.overlap))
#mcols(ltrs.dj)$ltrsim=splitAsList(mcols(ltrs.mask)$ltr_similarity[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$ltrsim=splitAsList(mcols(ltrs.mask)$ltr_similarity[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$ltrsim.c=sapply(mcols(ltrs.dj)$ltrsim, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$ltrsim.self=sapply(mcols(ltrs.dj)$ltrsim.c, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])

mcols(ltrs.dj)$geneorder=splitAsList(prot.mask$geneorder[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$geneorder.c=sapply(mcols(ltrs.dj)$geneorder, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$geneorder.self=sapply(mcols(ltrs.dj)$geneorder.c, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])


mcols(ltrs.dj)$nestlevel=splitAsList(mcols(ltrs.mask)$overlaps[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$nestlevel=sapply(mcols(ltrs.dj)$nestlevel, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$nestlevel=sapply(mcols(ltrs.dj)$nestlevel, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])




mcols(ltrs.dj)$ownname=splitAsList(mcols(ltrs.mask)$name[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$ownname=sapply(mcols(ltrs.dj)$ownname, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$ownname=sapply(mcols(ltrs.dj)$ownname, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])


mcols(ltrs.dj)$family=splitAsList(mcols(ltrs.mask)$family[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$family=sapply(mcols(ltrs.dj)$family, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$family=sapply(mcols(ltrs.dj)$family, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])





ltrs.overlap.equal=findOverlaps(ltrs.dj, ltrs.mask, type='equal')
ltrs.gff=data.frame(seqnames(ltrs.dj), '.', 'region', start(ltrs.dj), end(ltrs.dj), '.', '+', '.', 'col9')
names(ltrs.gff)=c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')
ltrs.gff$V3=as.character(ltrs.gff$V3)
ltrs.gff$V9=as.character(ltrs.gff$V9)
mcols(ltrs.dj)$value=sapply(mcols(ltrs.dj)$value, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$nestlevels=str_count(mcols(ltrs.dj)$value, ',')+1
for (linenum in 1:length(ltrs.dj)){
  if ( linenum %in% queryHits(ltrs.overlap.equal)){
    ltrs.gff[linenum,3]='LTR_retrotransposon'
    ltrs.gff[linenum,9]=paste('ID=', mcols(ltrs.dj)$ownname[linenum], 
#    					'; Parent=', mcols(ltrs.dj)$value[linenum], 
    					'; Name=NEST0,', mcols(ltrs.dj)$family[linenum],',', mcols(ltrs.dj)$ltrsim.self[linenum], ',', mcols(ltrs.dj)$geneorder.self[linenum], ',overlappingTEcopies', mcols(ltrs.dj)$nestlevel[linenum], ',TEsthisiswithin', mcols(ltrs.dj)$nestlevels[linenum],
    					'; Alias=', mcols(ltrs.dj)$ltrsim.c[linenum], ',', mcols(ltrs.dj)$geneorder.c[linenum], sep='')
  }
  else{
    ltrs.gff[linenum,3]='transposon_fragment'
    ltrs.gff[linenum,9]=paste('ID=FRAG_', mcols(ltrs.dj)$ownname[linenum], 
#    					'; Parent=', mcols(ltrs.dj)$value[linenum], 
    					'; Name=NEST0,', mcols(ltrs.dj)$family[linenum],',', mcols(ltrs.dj)$ltrsim.self[linenum], ',', mcols(ltrs.dj)$geneorder.self[linenum], ',overlappingTEcopies', mcols(ltrs.dj)$nestlevel[linenum], ',TEsthisiswithin', mcols(ltrs.dj)$nestlevels[linenum],
    					'; Alias=', mcols(ltrs.dj)$ltrsim.c[linenum], ',', mcols(ltrs.dj)$geneorder.c[linenum], sep='')
  }
}




write.table(ltrs.gff, paste(GENOMENAME, 'nested_subtracted_ltr.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)




