#!/bin/bash -login
#SBATCH -D /group/jrigrp8/mstitzer/multiple_maize_assemblies/w22/helitron/
#SBATCH -o /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/cluster_families-stdout-%j.txt
#SBATCH -e /group/jrigrp8/mstitzer/multiple_maize_assemblies/slurm-log/cluster_families-stderr-%j.txt
#SBATCH -J cluster_helitron
set -e
set -u


CPU=1

#USEARCH=~/software/usearch8.0.1623_i86linux32 
USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
HELFILE=W22.Helitron.noB73v4Fam.fa
HELBASE=$( basename $HELFILE .fa )

${USEARCH} -allpairs_global ${HELFILE} -blast6out ${HELBASE}.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

${SILIX} ${HELFILE} ${HELBASE}.allvall.8080.out -f HEL -i 0.8 -r 0.8 --net > ${HELBASE}.8080.fnodes


