#!/bin/bash -login
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr/families
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/add_fams-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/add_fams-stderr-%j.txt
#SBATCH -J add_fam
set -e
set -u

export PATH=$PATH:~/software/bin/

CPU=64

#USEARCH=~/software/usearch8.0.1623_i86linux32 
#USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
#SILIX=~/software/bin/silix
#LTRFILE=B73V4.both_pseudo_AND_unplaced.allnested.5ltr.newNamewFam.fa
#NEWTEFILE=W22.5ltr.fa
#NEWTEFILE=zea_diploperennis.5ltr.fa
#i=$SLURM_ARRAY_TASK_ID
#NEWTEFILE=B104.5ltr.fa

source ../CONFIG.sh

cat ../ltrdigest/${GENOME}.ltrdigest_5ltr.fas ../mask_subtract/subtract*_${GENOME}/ltrdigest/${GENOME}.*5ltr.fas > ${GENOME}.5ltr.fa


NEWTEFILE=${GENOME}.5ltr.fa

samtools faidx $NEWTEFILE

#mkdir -p subtracted
if [ ! -f ${NEWTEFILE}.MCSnames.id80.cov80.out ]
then
### this adds new TEs to existing B73 families
$USEARCH --usearch_global $NEWTEFILE -db $LTRFILE -id 0.8 -blast6out ${NEWTEFILE}.MCSnames.id80.cov80.out -strand both -query_cov 0.8 -target_cov 0.8 -top_hits_only --threads $CPU
fi

## the following are not generalized yet!
Rscript get_nofams.R ## this has hardcoded in to the r script w22 specific things!! 

xargs samtools faidx ${NEWTEFILE} < ${GENOME}.noExistingFam.txt >> ${GENOME}.noExistingFam.fa

${USEARCH} -allpairs_global ${GENOME}.noExistingFam.fa -blast6out ${GENOME}.noExistingFam.allvall.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads 1

${SILIX} ${GENOME}.noExistingFam.fa ${GENOME}.noExistingFam.allvall.out -f LTR -i 0.8 -r 0.8 --net > ${GENOME}.noExistingFam.8080.fnodes


