#!/bin/bash -login
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr/mask_subtract
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/run_nested_analysis-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/run_nested_analysis-stderr-%j.txt
#SBATCH -J ltrdigest_nest
set -e
set -u

source ../CONFIG.sh

GENOMETOOLS=$GENOMETOOLS151

export PATH=$PATH:${SOFTWAREPATH}

### genome tools path
#GENOMETOOLS=/home/mstitzer/software/genometools-1.5.1/bin/gt
GENOMETOOLS=/home/mstitzer/software/genometools-1.5.7/bin/gt
GENOMEBASE=W22__Ver12


i=$SLURM_ARRAY_TASK_ID
SUBTRACT=${GENOME}.subtract${i}

#GENOMEFASTA=../../${GENOME}.fasta

#MEMLIM=96GB
CPU=1



###################
## run ltrdigest ##
###################

mkdir -p ${OUTDIR}/mask_subtract/subtract${i}_${GENOME}/ltrdigest


echo $GENOME.ltrdigest
echo $GENOME.ltrharvest.sorted.gff3
echo $GENOME
echo $GENOME.ltrdigest.gff3

#$GENOMETOOLS -j $CPU ltrdigest -outfileprefix hardmask${i}/ltrdigest/$GENOME.ltrdigest -v yes -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- hardmask${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME > hardmask${i}/$GENOME.ltrdigest.gff3
## troubleshooting thinking multiple processors weren't working, but looks like it's storing stuff in /tmp on the compute node

### fix the .des files that look funny
#gt encseq encode -des yes -dna yes B73.Mhap2.quiver.subtract13.fa
#$GENOMETOOLS encseq encode -des yes -ssp no -sds no -md5 no -dna yes -indexname $GENOME $GENOMEFASTA


$GENOMETOOLS -j $CPU ltrdigest -outfileprefix ${OUTDIR}/mask_subtract/w22_mask_subtract/subtract${i}/ltrdigest/$SUBTRACT.ltrdigest -trnas ../eukaryotic-tRNAs.fa -hmms ../gydb_hmms/GyDB_collection/profiles/*.hmm -- ${OUTDIR}/mask_subtract/w22_mask_subtract/subtract${i}/$SUBTRACT.ltrharvest.sorted.gff3 $SUBTRACT > ${OUTDIR}/mask_subtract/w22_mask_subtract/subtract${i}/$SUBTRACT.ltrdigest.gff3

## troubleshooting further -- most basic command
##### I DO NOT KNOW WHY BUT I HAVE TO USE genometools-1.5.1 to get the hmm domains to work. I'm giving up on this for now. 

#$GENOMETOOLS -j $CPU ltrdigest -outfileprefix digesttest -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -trnas eukaryotic-tRNAs.fa hardmask${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME
#$GENOMETOOLS -j $CPU ltrdigest -outfileprefix digesttest -hmms all_gydb_profiles.hmmer3b.hmm -trnas eukaryotic-tRNAs.fa hardmask${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME

