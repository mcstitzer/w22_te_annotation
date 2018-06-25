#!/bin/bash -login
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/run_nested_analysis-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/run_nested_analysis-stderr-%j.txt
#SBATCH -J ltrharvest
#SBATCH -p bigmemh
#SBATCH --mem=96000
set -e
set -u

source CONFIG.sh

export PATH=$PATH:$SOFTWAREPATH

### genome tools path
GENOMETOOLS=$GTPATH

# name the file stem based on suffixator index
GENOME=W22__Ver12
GENOMEFASTA=../${GENOME}.genome.fasta

#MEMLIM=64GB
#CPU=8




mkdir -p $OUTDIR
####################################
## get protein domain hmms, tRNAs ##
####################################

if [ ! -d gydb_hmms ]
then
	./get_tRNA_hmm_dbs.sh
fi

###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################

if [ ! -f ${GENOME}.suf ]
then
$GENOMETOOLS suffixerator -db $GENOMEFASTA -indexname $GENOME -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
fi

#####################
## Run LTR harvest ##
#####################

mkdir -p ${OUTDIR}/outinner
if [ ! -f ${OUTDIR}/${GENOME}.ltrharvest.sorted.gff3 ]
then

## all defaults except for maxdistltr (default 15000)
$GENOMETOOLS ltrharvest -index $GENOME -gff3 ${OUTDIR}/$GENOME.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr 20000 -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner ${OUTDIR}/outinner/${GENOME}.ltrharvest.outinner.fa -out ${OUTDIR}/${GENOME}.ltrharvest.fa > ${OUTDIR}/${GENOME}.ltrharvest.out

$GENOMETOOLS gff3 -sort ${OUTDIR}/$GENOME.ltrharvest.gff3 > ${OUTDIR}/$GENOME.ltrharvest.sorted.gff3

fi

###################
## run ltrdigest ##
###################

mkdir -p ${OUTDIR}/ltrdigest

if [ ! -f ${OUTDIR}/${GENOME}.ltrdigest.gff3 ]
then
$GENOMETOOLS -j $CPU ltrdigest -outfileprefix ${OUTDIR}/ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- ${OUTDIR}/$GENOME.ltrharvest.sorted.gff3 $GENOME > ${OUTDIR}/$GENOME.ltrdigest.gff3
fi

