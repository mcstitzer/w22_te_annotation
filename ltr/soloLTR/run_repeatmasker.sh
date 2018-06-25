#!/bin/bash -login
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr/soloLTR
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/repeatmasker-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/repeatmasker-stderr-%j.txt
#SBATCH -J repeatmasker
set -e
set -u


export PATH=$PATH:/home/mstitzer/software/bin

## repeatmasker will make output files in whatever directory the orignial genome is in
## copy it over to keep them here.
GENOME=W22__Ver12.$1
#GENOMEFASTA=../../${GENOME}.fasta
GENOMEFASTA=${GENOME}.fasta

#cp $GENOMEFASTA .
#LTRFILE=../families/B73V4.both_pseudo_AND_unplaced.allnested.5ltr.newNamewFam.fa

LTRFILE=b73_w22_5ltr.MCSname.fa
if [ ! -f $LTRFILE ]
then
cat ../families/B73V4.both_pseudo_AND_unplaced.allnested.5ltr.newNamewFam.fa ../../w22/ltr/W22.5ltr.MCSnames.fa > b73_w22_5ltr.MCSname.fa
fi


#python switch_fasta_names.py ../families/${GENOME}.allnested_5ltr.fa ../arabidopsis_OrigName_TAIR10.txt > ${GENOME}.allnested_5ltr.newNamewFam.fa
#LTRFILE=${GENOME}.allnested_5ltr.newNamewFam.fa
#TAIR10_Chr.all.allnested_5ltr.fa



module load trf crossmatch rmblast hmmer repeatmasker


#~/software/RepeatMasker/RepeatMasker -pa 32 -lib $LTRFILE -qq -gff -nolow -no_is ${GENOME}.fasta
#RepeatMasker -pa 24 -lib $LTRFILE -qq -gff -nolow -no_is ${GENOME}.fasta
~/software/RepeatMasker_4-0-7/RepeatMasker -pa 24 -lib $LTRFILE -qq -gff -nolow -no_is ${GENOME}.fasta


Rscript repeatmasker_gff_to_flanking_gff.R


bedtools getfasta -fi ${GENOME}.fasta -bed ${GENOME}.fasta.up.gff3 -tab -fo ${GENOME}.fasta.up.tabout
bedtools getfasta -fi ${GENOME}.fasta -bed ${GENOME}.fasta.down.gff3 -tab -fo ${GENOME}.fasta.down.tabout
