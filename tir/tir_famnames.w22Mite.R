library(rtracklayer)
library(stringr)

tir=read.table('w22_alltir.plusdetectMITE.tab')

names(tir)=c('mtec', 'chr', 'start', 'end', 'tirscore', 'V6', 'V7', 'TIR1', 'TIR2')
tir.gr=GRanges(seqnames=tir$chr, ranges=IRanges(start=tir$start, end=tir$end))
mcols(tir.gr)$score=tir$tirscore
selfOver=findOverlaps(tir.gr, drop.self=T, ignore.strand=T)
rmRows=sapply(1:length(selfOver), function(x){ 
scoreA=mcols(tir.gr)$score[queryHits(selfOver)[x]]
scoreB=mcols(tir.gr)$score[subjectHits(selfOver)[x]]
if(scoreA>=scoreB){   ## arbitrarily prioritize the first entry if the scores are equal.
return(queryHits(selfOver)[x])
}else if (scoreB>scoreA){
return(subjectHits(selfOver)[x])
}
}
)


tir=tir[-rmRows,]
tir.gr=tir.gr[-rmRows,]

tir.gr$mtec=tir$mtec
tir.gr$mtecfamnum=substr(tir.gr$mtec, 7,11)
tir.gr$famname=paste(tir.gr$superfam, tir.gr$mtecfamnum, sep='')

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
mcols(tir.gr)$Name=NULL
for (x in names(table(mcols(tir.gr)$famname))){
  mcols(tir.gr)$Name[mcols(tir.gr)$famname==x & !is.na(mcols(tir.gr)$famname)]=paste('B73v4', str_pad(1:sum(mcols(tir.gr)$famname[!is.na(mcols(tir.gr)$famname)]==x), 5, pad='0'),sep='')
}

tir.gr$sup=substr(tir.gr$mtec, 1,3)
tir.gr$ID=paste0(tir.gr$sup, tir.gr$Name)

d=data.frame(chr=seqnames(tir.gr), source='TARGeT', kind='terminal_inverted_repeat_element', start(tir.gr), end(tir.gr), '.', '*', '.', Name=tir.gr$ID)
write.table(d, 'w22_tir_copies.withdetectMITE.gff3', sep='\t', row.names=F, col.names=F, quote=F)

