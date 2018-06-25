#!/bin/bash
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/w22/tir/target_mTEA
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/target_dna-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/target_dna-stderr-%j.txt
#SBATCH -J target_mite
#SBATCH -p med


#NO SBATCH --ntasks-per-node 16
#NO SBATCH --mem 24000

GENOME=../../../W22__Ver12.genome.fasta

module load blast/2.2.26
export PATH=$PATH:/home/mstitzer/software/bin
export PATH=$PATH:/home/mstitzer/projects/agpv4_te_annotation/tir/mTEA/lib/blogo/
#export $PERL5LIB
export PERL5LIB=$PERL5LIB:/home/mstitzer/projects/agpv4_te_annotation/tir/mTEA/lib/blogo/

#FILES=($(ls -1 ~/te_reference_fasta/individual_fasta/DTC_*))

readarray -t MITEs < detectMITE_positions_tirtsd.txt

#IFS='\t'

#MITE=${FILES[$SLURM_ARRAY_TASK_ID]}
#echo $MITE

lineNumber=$SLURM_ARRAY_TASK_ID+1
echo $lineNumber
#read -r -a MITE <<< "${MITEs[$SLURM_ARRAY_TASK_ID]}"
#MITE=$(awk 'NR==$lineNumber' detectMITE_positions_tirtsd.txt) # | while read line; do

MITE=($(echo ${MITEs[$lineNumber]}))
echo $MITE
DNATE=${MITE[3]}
TSD=${MITE[4]}
TIR=${MITE[5]}
#done
echo $MITE
echo $DNATE
echo $TIR
echo $TSD

### set up a file name for each individual TE fasta
#FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}
#echo $FILENAME

### RUN TARGeT
#FILE=$(basename "$FILENAME")
#echo $FILE
#DNATE="${FILE%.*}"
#echo $DNATE
#DNATE=${MITE[4]}
#TSD=${MITE[5]}
#TIR=${MITE[6]}


### when run out of memory, use this with more memory next time?
#if ${LTR}.flank.fa.tab not in ls:

if [ ! -f ${DNATE}.tir.fa.tab ]
then

mkdir -p $DNATE
python ~/software/TARGeT/target.py -q ../b73v4_tir/${DNATE}.fa -t nucl -o $DNATE -i s -P 1 -S PHI -DB -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.3 $GENOME ${DNATE}_target
### db generation for the first time to index the genome
#python ~/software/TARGeT/target.py -q $FILENAME -t nucl -o $DNATE -i s -P 1 -S PHI -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.3 $GENOME ${DNATE}_target


### Convert names of flanking fasta file
if [ ! -f ${DNATE}.flank_adj ]
then
python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank > ${DNATE}.flank.fa
else
python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank_adj > ${DNATE}.flank.fa
fi



### RUN mTEA

sup=${DNATE:0:3}
echo $sup

if [ "$sup" == "DTC" ]
then
## CACTA : -c NNN -t CACTNNNNNNNNN -d 3 -s 0
perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNN -t CACTNNNNNNNNN -d 3 -s 0
elif [ "$sup" == "DTA" ]
then
## hAT : -c NNNNNNNN -t NNNNNNNNNNN -d 2 -s 1
perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNNNNNNN -t NNNNNNNNNNN -d 2 -s 1
elif [ "$sup" == "$DTM" ]
then
## Mutator : -c NNNNNNNNNNN -t NNNNNNNNNNN -d 4 -s 1
perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNN -t CACTNNNNNNNNN -d 3 -s 0
elif [ "$sup" == "DTT" ]
then
## Tc1/Mariner : -c TA -t NNNNNNNNNNNN -d 2 -s 0
perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c TA -t NNNNNNNNNNNN -d 2 -s 0
elif [ "$sup" == "DTH" ]
then
## PIF/Harbinger : -c TNN -t NNNNNNNNNNNNNN -d 3 -s 0
perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c TNN -t NNNNNNNNNNNNNN -d 3 -s 0
else ## this means it is a DTX, so we need to get TSD and TIR from the file
# here, we're using TSD and TIR lengths, and generating a string of N's the appropriate length
# and using mismatches by rounding to get 80% identity between TSD and TIR
perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c $(printf -v sep '%*s' $TSD ; echo "${sep// /N}") -t $(printf -v sep '%*s' $TIR ; echo "${sep// /N}") -d $(($TIR*0.2)) -s $(($TSD*0.2))

fi

#rm -rf $DNATE
fi

### so to submit, grab files and export them
####### DO ALL OF THIS IN THE SHELL YOU'RE ABOUT TO SUBMIT IN
### export FILES=($(ls -1 DTA*))
### NUMFILES=${#FILES[@]}
### ARRAYNUM=$(($NUMFILES - 1)) ## because slurm array is zero based

## if [ $ARRAYNUM -ge 0 ]; then
## sbatch --array=0-$ARRAYNUM tir_search_genome.sh
## fi
