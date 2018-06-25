#!/bin/bash -login

GENOME=W22
GENOMEFASTA=../W22__Ver12.genome.fasta
CPU=1

## path to SILIX for clustering
SILIX=/home/mstitzer/software/bin/silix
## path to VSEARCH for matching sequences
VSEARCH=/home/mstitzer/software/vsearch/bin/vsearch

#### get the SINE-Finder program from Wenke et al. 2011
wget http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt 

#### change name
mv Supplemental_Data_Set_1-sine_finder.txt sine_finder.py

#### run sinefinder
#### I haven't been able to get sine_finder to work with reverse sequences, as it seems to report TSDs wrong on the reverse strand.
####   so I'm only reporting on the forward strand.
### -f both : outputs csv and fasta
python sine_finder.py -T chunkwise -V1 -f both -o F ../${GENOMEFASTA}

#### sine_finder outputs the fasta with the TSD included. I remove these here, so they aren't considered when clustering into families
mv ../${GENOME}-matches.fasta .
mv ../${GENOME}-matches.csv .
python remove_tsd_sinefinder.py ${GENOME}-matches.fasta ${GENOME}-matches.noTSD.fa

#### vsearch to identify homology, silix to cluster
$VSEARCH --usearch_global ${GENOME}-matches.noTSD.fa -db RST.B73v4.fa -id 0.8 -blast6out ${GENOME}-matches.noTSD.MCSnames.id80.cov80.out -query_cov 0.8 -target_cov 0.8 --threads $CPU --minseqlength 1


## the following are not generalized yet!
samtools faidx ${GENOME}-matches.noTSD.fa
Rscript get_nofams.R ## this has hardcoded in to the r script w22 specific things!! 

xargs samtools faidx ${GENOME}-matches.noTSD.fa < W22.RST.noB73v4Fam.txt >> W22.RST.noB73v4Fam.fa

${VSEARCH} -allpairs_global W22.RST.noB73v4Fam.fa -blast6out W22.RST.noB73v4Fam.allvall.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads 1


# single linkage cluster those that are 80% identical to each other.
$SILIX W22.RST.noB73v4Fam.fa W22.RST.noB73v4Fam.allvall.out -f SINE -i 0.8 -r 0.8 > ${GENOME}.RST.noB73v4fam.8080.fnodes

### cluster into families and output final gff with this R script
Rscript generate_gff_SINE.R $GENOME


