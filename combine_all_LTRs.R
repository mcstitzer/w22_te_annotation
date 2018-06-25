library(rtracklayer)
library(stringr)
library(plyr)



setwd('../nested_ltr/')

GENOMENAME='W22__Ver12'
SHORTID='Zm00004b'
MAXSUBTRACT=80

## use the sorted one here to not run into problems with ordering not matching digest results
masked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '_', GENOMENAME, '/', GENOMENAME, '.subtract', 1:MAXSUBTRACT, '.ltrharvest.contignames.sorted.gff3.contigpositions.sorted.gff3', sep='')
filenames=c(paste('mask_subtract/', GENOMENAME, '.ltrharvest.contignames.sorted.gff3', sep=''), masked)

#filenames=filenames[-c(30,59)]

#### read in all the masked TEs


#ltrs=lapply(filenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
ltrs=lapply(filenames, function(x) {a=import.gff(x) ## change in rtracklayer, version='3')  ### no asRangedData in new version of R
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
##### MAKE SURE THESE HAVE THE ^H REMOVED!

dignames=paste('families/',GENOMENAME, '_tabout/', GENOMENAME, '.subtract', 1:MAXSUBTRACT, '.ltrdigest_tabout.csv', sep='')
digfiles=c(paste('ltrdigest/', GENOMENAME, '.ltrdigest_tabout.csv', sep=''), dignames)

#digfiles=digfiles[-c(30,59)]


split.fam=function(s) strsplit(as.character(s), '_')[[1]][1]
prot=lapply(digfiles, function(x) {m=read.table(x, sep='\t', header=T)
                                                                        m$sequence=str_split_fixed(as.character(m$sequence), '_', 3)[,1]
                                                                        m$rle=''
                                                                        print(x)
                                                                        if (sum(m$Protein.domain.hits!='', na.rm=T)>0){
                                                                                m$rle[m$Protein.domain.hits!='']=sapply(which(m$Protein.domain.hits!=''), function(x)  rle(sapply(ldply(strsplit(as.character(m$Protein.domain.hits[x]), '/')), split.fam))$values)
                                                                                m$geneorder=sapply(1:nrow(m), function(x) paste(m$rle[[x]], collapse='.'))
                                                                                m$nestlevel=which(digfiles==x)-1
                                                                                } else{
                                                                                m$geneorder=''
                                                                                m$nestlevel=which(digfiles==x)-1
                                                                                }
                                                                        return(m)
                                                                        }
                        )

prot.mask=do.call(rbind, prot)


############# get the same for the ltrdigest gff3's so that I can add strand
dmasked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '_', GENOMENAME, '/', GENOMENAME, '.subtract', 1:MAXSUBTRACT, '.ltrdigest.gff3', sep='')
dfilenames=c(paste(GENOMENAME, '.ltrdigest.gff3', sep=''), dmasked)

#dfilenames=dfilenames[-c(30,59)]															
															
dltrs=lapply(dfilenames, function(x) {a=import.gff3(x) ## change in rtracklayer #, version='3', asRangedData=F)
									print(x)
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    ltr_similarity=as.numeric(mcols(a)$ltr_similarity) ## can use this to confirm same order
                                    nestlevel=which(dfilenames==x)-1
                                    ID=mcols(a)$ID
                                    mcols(a)=NULL ### get rid of this because what comes first determines columns, which differ between file
                                    mcols(a)$ltr_similarity=ltr_similarity
                                    mcols(a)$nestlevel=nestlevel
                                    mcols(a)$ID=ID
                                    return (a)
                                    }
            )
            
dltrs.mask=do.call(c, dltrs)


########## OH NO THESE ARE THE WRONG ORDER!!!!!######
###### these will be good to switch the order, as the start and end positions are the subtracted coordinates

### get ltrharvest results - need to do chromosomes and unplaced separately because i did them separately (and won't again!)

smasked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '_', GENOMENAME, '/', GENOMENAME, '.subtract', 1:MAXSUBTRACT, '.ltrharvest.gff3', sep='')
sfilenames=c(paste(GENOMENAME, '.ltrharvest.gff3', sep=''), smasked)

#sfilenames=sfilenames[-c(30,59)]															

sltrs=lapply(sfilenames, function(x) {a=import.gff3(x) # change in rtracklayer, version='3', asRangedData=F)
									print(x)
                                    mcols(a)$nest=which(sfilenames==x)-1
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    mcols(a)$ltr_similarity=as.numeric(mcols(a)$ltr_similarity)
                                    mcols(a)$nestlevel=which(sfilenames==x)-1
                                    return (a)
                                    }
            )
            
sltrs.mask=do.call(c, sltrs)


#########

## confirm these are in the same order

all(mcols(ltrs.mask)$ltr_similarity==mcols(dltrs.mask)$ltr_similarity)   ## concerning that this is not true
## weird behavior of genome tools in that it will reorder the ID number so each gff goes 1 to 
all(mcols(ltrs.mask)$ID==mcols(dltrs.mask)$ID)   ### but this has to be right! does this mean i'm switching contig names wrong?

## confirm same number of TEs in ltrs.mask and prots.mask ((note sltrs.mask and dltrs.mask might be different if diff seqname was at bottom of fasta))															
all(sort(table(seqnames(ltrs.mask)))==sort(table(prot.mask$sequence)))															
															
															
##fine now
### there are different orders in the digest and harvest and i don't know why but i'm solving it.
#mcols(ltrs.mask)$compare=paste('seq', mcols(ltrs.mask)$seq_number, mcols(ltrs.mask)$ltr_similarity, mcols(ltrs.mask)$nestlevel, sep='')
#mcols(dltrs.mask)$compare=paste(seqnames(dltrs.mask), mcols(dltrs.mask)$ltr_similarity, mcols(dltrs.mask)$nestlevel, sep='')
#ltrs.mask=ltrs.mask[order(match(mcols(ltrs.mask)$compare, mcols(dltrs.mask)$compare)),]




## these should be the same order as the ltrs.mask object
## also check that the start values are the same to be super sure
length(ltrs.mask)==nrow(prot.mask)
all(start(ltrs.mask)==prot.mask$element.start) ### this will be false now that these are subtracted!!!!
## if in different orders, this will be false (not sorted)
all(as.character(seqnames(ltrs.mask))==prot.mask$sequence)
all(as.character(seqnames(ltrs.mask))==prot.mask$sequence)


## if these are all fine, it's time to assign strand to the ltrs!

strand(ltrs.mask)=strand(dltrs.mask)



##############################################
##  Read in and process family assignments  ##
##############################################

######### Read in silix results
a=read.table('families/W22.5ltr.fa.MCSnames.id80.cov80.out')
a$V1=as.character(a$V1)
a$V2=as.character(a$V2)
a$fam=substr(a$V2, 1,8)


### get switch positions for chromosome to subtracted chromosome here

glen=readRDS(paste('mask_subtract/', GENOMENAME, '_gl_', MAXSUBTRACT, '.RDS', sep=''))



															
###########################################
##   Assign families and superfamilies   ##
###########################################

															
###  SUPER SUPER concerning that the ltrs.mask seqnames are wrong!! I switched to using the prot.mask names, which are appropriately switched.
															
mcols(ltrs.mask)$name=paste(seqnames(ltrs.mask), start(ltrs.mask), end(ltrs.mask), sep='_')
### this is the original name with positions on the subtracted chromosomes (not the coordinates from the full v4 chromosomes)
prot.mask$origname=paste(prot.mask$sequence, prot.mask$element.start, prot.mask$element.end, sep='_')
#prot.mask$origname[ltrs.mask$nestlevel>0]=paste(seqnames(ltrs.mask), prot.mask$element.start+1, prot.mask$element.end, sep='_')[ltrs.mask$nestlevel>0]
mcols(ltrs.mask)$origname=prot.mask$origname

## add the family column - different for the original intact level and the subtracted ones
mcols(ltrs.mask)$family=NA
mcols(ltrs.mask)$family=mapvalues(mcols(ltrs.mask)$origname, from=a$V1, to=a$fam, warn_missing=F)
ltrs.mask$family[substr(ltrs.mask$family,1,1)!='R']=NA
#mcols(ltrs.mask)$family[mcols(ltrs.mask)$nestlevel==0]=mapvalues(mcols(ltrs.mask)$name[mcols(ltrs.mask)$nestlevel==0], from=a$V2, to=a$V1)
#mcols(ltrs.mask)$family[mcols(ltrs.mask)$name %in% allfam$V1]=mapvalues(mcols(ltrs.mask)$name[mcols(ltrs.mask)$name %in% allfam$V1], from=allfam$V1, to=allfam$V2)


															
## for each copy, assign whether it can be attributed to either superfamily based on protein order
copia=c('GAG.AP.INT.RT.RNaseH.ENV', 'GAG.AP.INT.RT.RNaseH', 'INT.RT.RNaseH', 'AP.INT.RT.RNaseH', 'GAG.GAGCOAT.GAG.AP.INT.RT.RNaseH.ENV', 'GAG.AP.RT.RNaseH', 'AP.INT.RT.RNaseH.ENV', 'GAG.GAGCOAT.GAG.AP.INT.RT.RNaseH', 'INT.RT
.RNaseH.ENV')
gypsy=c('GAG.AP.RT.RNaseH.INT', 'GAG.AP.RT.RNaseH.INT.CHR', 'GAG.GAGCOAT.AP.RT.RNaseH.INT.CHR', 'RT.RNaseH.INT.CHR', 'GAG.GAGCOAT.GAG.AP.RT.RNaseH.INT.CHR', 'RT.RNaseH.INT', 'GAG.AP.RT.INT', 'AP.RT.RNaseH.INT.CHR', 'GAG.GAGCOAT.AP.RT.RNaseH.INT', 'GAG.RT.RNaseH.INT.CHR', 'RNaseH.INT', 'GAG.RT.RNaseH.INT', 'GAG.AP.RT.INT.CHR', 'RNaseH.INT.CHR', 'AP.RT.RNaseH.INT')
mcols(ltrs.mask)$superfam=NA
mcols(ltrs.mask)$superfam[prot.mask$geneorder %in% copia]='RLC'
mcols(ltrs.mask)$superfam[prot.mask$geneorder %in% gypsy]='RLG'
															
													
###########################
## Read new W22 families ##
###########################
															
															
															
wfam=read.table('../nested_ltr/families/W22.noB73v4Fam.8080.fnodes', stringsAsFactors=F)
ltrs.mask$w22fam=mapvalues(ltrs.mask$origname, from=wfam$V2, to=wfam$V1)
ltrs.mask$w22fam[substr(ltrs.mask$w22fam, 1,1)!='L']=NA
															
## assign superfamily to each family based on whether protein order can be assigned
fam.sup=sapply(names(table(mcols(ltrs.mask)$w22fam)), function(x) names(which.max(table(mcols(ltrs.mask)$superfam[which(mcols(ltrs.mask)$w22fam==x)]))))
mcols(ltrs.mask)$w22sup=mapvalues(mcols(ltrs.mask)$w22fam, from=names(unlist(fam.sup)), to=unlist(fam.sup))
mcols(ltrs.mask)$w22sup[!mcols(ltrs.mask)$w22sup %in% c('RLC', 'RLG') & !is.na(ltrs.mask$w22fam)]='RLX'  ## only do for those that already have a family!!

ltrs.mask$w22sup[is.na(ltrs.mask$w22fam)]=NA
	       
### switch the LTR naming of families from silix to one that reflects superfamily ala Wicker et al. (2007), and is 5 digits long
#mcols(ltrs.mask)$famname=sapply(1:length(ltrs.mask), function(x) sub('LTR0', mcols(ltrs.mask)$fam.sup[x], mcols(ltrs.mask)$family[x]))
															
															
#############################################
## To update existing RLX to a superfamily ##
#############################################
				
##not done yet]
ltrs.mask$b73sup=substr(ltrs.mask$family,1,3)
sum(ltrs.mask$b73sup=='RLX' & ltrs.mask$superfam %in% c('RLC', 'RLG'), na.rm=T)
tail(sort(table(ltrs.mask$family[ltrs.mask$b73sup=='RLX' & ltrs.mask$superfam %in% c('RLC', 'RLG')])))
length(sort(table(ltrs.mask$family[ltrs.mask$b73sup=='RLX' & ltrs.mask$superfam %in% c('RLC', 'RLG')])))
## so there are 626 of these!
	
	      
############################
### Assign families      ###
############################	
## get all possible b73 family names!
##
b=import.gff3('~/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73v4.TE.gff3')
b$fam=substr(b$ID, 1,8)
b$sup=substr(b$ID, 1,3)
maxFamName=max(as.numeric(substr(b$fam,4,8)[b$sup %in% c('RLG', 'RLX', 'RLC')]))	       

## This gave the same answer but could have been terrible if the largest numbered fam wasn't in W22	       
#maxFamName=max(as.numeric(substr(ltrs.mask$family, 4,8)), na.rm=T)
	       
### switch the LTR naming of families from silix to one that reflects superfamily ala Wicker et al. (2007), and is 5 digits long
#mcols(ltrs.mask)$famname=sapply(1:length(ltrs.mask), function(x) sub('LTR0', mcols(ltrs.mask)$w22sup[x], mcols(ltrs.mask)$w22fam[x]))

## first for new w22 families  ## this is super slow!!
for (x in 1:length(table(mcols(ltrs.mask)$w22fam))){
  famNum=maxFamName + x  # we want to increment families
  famName=names(rev(sort(table(mcols(ltrs.mask)$w22fam))))[x]  ## and keep track of original families
  mcols(ltrs.mask)$famname[mcols(ltrs.mask)$w22fam==famName & !is.na(mcols(ltrs.mask)$w22fam)]=paste(famNum, SHORTID, str_pad(1:sum(mcols(ltrs.mask)$w22fam[!is.na(mcols(ltrs.mask)$w22fam)]==famName), 5, pad='0'), sep='')
}	       
## assign rl superfam
mcols(ltrs.mask)$famnameW22=sapply(1:length(ltrs.mask), function(x) sub('^', mcols(ltrs.mask)$w22sup[x], mcols(ltrs.mask)$famname[x]))
	       
## assign 11 digit copy name (SHORTIDXXXXX) for each copy in an existing B73 family
ltrs.mask$Name=NA
for (x in names(table(mcols(ltrs.mask)$family))){
  mcols(ltrs.mask)$Name[mcols(ltrs.mask)$family==x & !is.na(mcols(ltrs.mask)$family)]=paste(x, SHORTID, str_pad(1:sum(mcols(ltrs.mask)$family[!is.na(mcols(ltrs.mask)$family)]==x), 5, pad='0'), sep='')
}
ltrs.mask$Name[is.na(ltrs.mask$b73sup)]=ltrs.mask$famnameW22[is.na(ltrs.mask$b73sup)]			

## done!!!!
sum(is.na(ltrs.mask$Name))
				
############################
### Keep EVERYTHING      ###
############################				
#### output a table so that I can recluster these. will need names for fasta headers and subtracted level
### write them separately for unplaced and chromosomes so that I can split them easily.
write.table(data.frame(prot.mask$nestlevel[is.na(mcols(ltrs.mask)$family) & prot.mask$sequence %in% as.character(1:10)], prot.mask$origname[is.na(mcols(ltrs.mask)$family)& prot.mask$sequence %in% as.character(1:10)]), paste0(GENOMENAME, '_copies_without_family_chromosome.txt'), quote=F, row.names=F, col.names=F, sep='\t')
write.table(data.frame(prot.mask$nestlevel[is.na(mcols(ltrs.mask)$family) & !prot.mask$sequence %in% as.character(1:10)], prot.mask$origname[is.na(mcols(ltrs.mask)$family)& !prot.mask$sequence %in% as.character(1:10)]), paste0(GENOMENAME, '_copies_without_family_unplaced.txt'), quote=F, row.names=F, col.names=F, sep='\t')
## this is for regular, all in one
##write.table(data.frame(prot.mask$nestlevel[is.na(mcols(ltrs.mask)$family)], prot.mask$origname[is.na(mcols(ltrs.mask)$family)]), 'copies_without_family.txt', quote=F, row.names=F, col.names=F, sep='\t')


### save these to update once i have families for all!
saveRDS(prot.mask, paste(GENOMENAME, '_prot.mask.dataframe.RDS', sep=''))
saveRDS(ltrs.mask, paste(GENOMENAME, '_ltrs.mask.genomicranges.RDS', sep=''))

				
				
## final gff

ltrs.gff=data.frame(seqnames(ltrs.mask), 'LTRharvest', 'LTR_retrotransposon', start(ltrs.mask), end(ltrs.mask), '.', strand(ltrs.mask), '.', paste('ID=', mcols(ltrs.mask)$Name, ';Name=', paste(mcols(ltrs.mask)$Name, paste('LTRsimilarity', mcols(ltrs.mask)$ltr_similarity, sep=''), sep='_'), sep=''))
write.table(ltrs.gff, paste(GENOMENAME, '.allFamilies.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

write.table(ltrs.gff[width(ltrs.mask)<=100000,], paste(GENOMENAME, '.allFamilies.100kbFilter.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

#saveRDS(ltrs.mask, paste(GENOMENAME, '.ltrs.mask.GRanges.RDS', sep=''))

## output original names and new names for easy comparison
write.table(data.frame(w22_orig=ltrs.mask$origname, MCSname=ltrs.mask$Name), 'original_name_to_MCSname.w22.txt', quote=F, sep='\t', row.names=F, col.names=T)
				
## output all families, with whether new or not
## need to update to use all B73 families!!!	
allFam=data.frame(fam=sort(unique(c(b$fam[b$sup %in% c('RLC', 'RLG', 'RLX')], substr(ltrs.mask$Name, 1,8)))), genome='b73')
allFam$genome=as.character(allFam$genome)
allFam[substr(allFam$fam,4,8)>maxFamName, 'genome']='w22'

write.table(allFam, 'all_families_added.w22.txt', quote=F, col.names=T, row.names=F, sep='\t')
				   
				   
				   
				   
				   
