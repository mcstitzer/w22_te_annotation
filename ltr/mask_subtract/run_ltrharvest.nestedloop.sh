#!/bin/bash -login
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr/mask_subtract
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/run_nested_analysis-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/run_nested_analysis-stderr-%j.txt
#SBATCH -J ltrharvest_nest
#SBATCH -p bigmemh
#SBATCH --mem=96000
set -e
set -u

source ../CONFIG.sh

mkdir -p $OUTDIR/mask_subtract

if [ ! -e $OUTDIR/mask_subtract/collapse_chromosomes.py ]
then
	cp collapse_chromosomes.py $OUTDIR/mask_subtract
fi

if [ ! -e $OUTDIR/mask_subtract/convert_ltrharvest_seq_gff_to_contignames.py ]
then
	cp convert_ltrharvest_seq_gff_to_contignames.py $OUTDIR/mask_subtract
fi

cd $OUTDIR/mask_subtract

export PATH=$PATH:${SOFTWAREPATH}

#MEMLIM=96GB

### genome tools path
#GENOMETOOLS=/home/mstitzer/software/genometools-1.5.7/bin/gt

i=1
#GENOME=W22__Ver12
#GENOME=PH207
#GENOME=B104
#GENOME=zea_diploperennis


## the first one is different, becuase need to set up subtract directory structure.

## this needs to be made more general so it is relative to OUTDIR
python convert_ltrharvest_seq_gff_to_contignames.py ../${GENOME}.ltrharvest.gff3 ../${GENOME}.des > ${GENOME}.ltrharvest.contignames.gff3
grep "LTR_retrotransposon	" ${GENOME}.ltrharvest.contignames.gff3 > ${GENOME}.ltrharvest.contignames.ltrretrotransposon.gff3
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" ${GENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > ${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx ${GENOMEFASTA}
bedtools complement -i ${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ${GENOMEFASTA}.fai > ${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
bedtools getfasta -fi ${GENOMEFASTA} -bed ${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo ${GENOME}.subtract1.fa

## concatenate the entries by chromosome
python collapse_chromosomes.py ${GENOME}.subtract1.fa > ${GENOME}.temp
mv ${GENOME}.temp ${GENOME}.subtract1.fa

### index this fasta
$GENOMETOOLS suffixerator -db ${GENOME}.subtract1.fa -indexname ${GENOME}.subtract1 -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
mkdir -p subtract1_${GENOME}
$GENOMETOOLS ltrharvest -index ${GENOME}.subtract1 -gff3 subtract1_${GENOME}/${GENOME}.subtract1.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr 21000 -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner subtract1_${GENOME}/${GENOME}.subtract1.ltrharvest.outinner.fa -out subtract1_${GENOME}/${GENOME}.subtract1.ltrharvest.fa > subtract1_${GENOME}/${GENOME}.subtract1.ltrharvest.out
gt gff3 -sort subtract1_${GENOME}/${GENOME}.subtract1.ltrharvest.gff3 > subtract1_${GENOME}/${GENOME}.subtract1.ltrharvest.sorted.gff3

### go until there are no more LTR TEs in this nested form 
#while [ grep -c ltr_retrotransposon ${GENOME}.hardmask${i}.ltrharvest.gff3 -gt 0 ]
while [ $i -le 80 ]
do

# name the file stem based on suffixator index

OLDINDEX=$i
OLDGENOME=${GENOME}.subtract${i}
i=$(( $i + 1 ))
NEWINDEX=$i
NEWGENOME=${GENOME}.subtract${i}

OLDGENOMEFASTA=${OLDGENOME}.fa
NEWGENOMEFASTA=${NEWGENOME}.fa


###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################


## switch genome tools back to their real contig names
python convert_ltrharvest_seq_gff_to_contignames.py subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.gff3 ${OLDGENOME}.des > subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.gff3
## only get the LTR_retrotransposon records to keep Rscript from repeating computation
grep "LTR_retrotransposon	" subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.gff3 > subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.ltrretrotransposon.gff3
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx ${OLDGENOMEFASTA}

###################UNIMPLEMENTED IDEA, this is now done after the fact!
### run the rscript: read in gff, read in RDS of genomelist, update genomelist with updatePos, write gff with changed positions


## find regions not covered by TEs
bedtools complement -i subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ${OLDGENOME}.fa.fai > subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
bedtools getfasta -fi $OLDGENOMEFASTA -bed subtract${OLDINDEX}_${GENOME}/${OLDGENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo $NEWGENOMEFASTA

## concatenate the entries by chromosome
python collapse_chromosomes.py $NEWGENOMEFASTA > ${NEWGENOMEFASTA}.temp
mv ${NEWGENOMEFASTA}.temp $NEWGENOMEFASTA

### index this fasta
$GENOMETOOLS suffixerator -db ${NEWGENOMEFASTA} -indexname ${NEWGENOME} -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM


#####################
## Run LTR harvest ##
#####################

mkdir -p subtract${NEWINDEX}_${GENOME}
### allow extra 1kb for each iteration, because we miss insertions that are not structural
MAXLEN=$(($i * 1000 + 20000))    ### so 20kb for the first hardmask, plus the additional 1kb per round
## all defaults except for maxdistltr (default 15000)
$GENOMETOOLS ltrharvest -index ${NEWGENOME} -gff3 subtract${NEWINDEX}_${GENOME}/${NEWGENOME}.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr $MAXLEN -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner subtract${NEWINDEX}_${GENOME}/${NEWGENOME}.ltrharvest.outinner.fa -out subtract${NEWINDEX}_${GENOME}/${NEWGENOME}.ltrharvest.fa > subtract${NEWINDEX}_${GENOME}/${NEWGENOME}.ltrharvest.out

$GENOMETOOLS gff3 -sort subtract${NEWINDEX}_${GENOME}/${NEWGENOME}.ltrharvest.gff3 > subtract${NEWINDEX}_${GENOME}/${NEWGENOME}.ltrharvest.sorted.gff3


done


####################### THIS IS NOT HOW I DO IT - instead, run each ltrdigest on every subtracted entry as an array

## can then run all the ltrdigest in array form on all the gffs
###################
## run ltrdigest ##
###################

#mkdir -p hardmask5/ltrdigest

#$GENOMETOOLS -j 16 ltrdigest -outfileprefix hardmask5/ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- hardmask5/$GENOME.ltrharvest.sorted.gff3 $GENOME > hardmask5/$GENOME.ltrdigest.gff3

