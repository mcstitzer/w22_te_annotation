#!/usr/bin/env bash

## probably a good idea to use entire path on these guys

GENOME=W22__Ver12
GENOMEFASTA=/group/jrigrp8/mstitzer/multiple_maize_assemblies/W22__Ver12.fasta
GENOMETOOLS=/home/mstitzer/software/genometools-1.5.7/bin/gt
GENOMETOOLS151=/home/mstitzer/software/genometools-1.5.1/bin/gt
OUTDIR=/group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr/ ## will make this directory if needed, nice to have outside of this github repository
## general directory of various required software noted in README.md (if not loading modules, etc)
SOFTWAREPATH=/home/mstitzer/software/bin

MEMLIM=64GB
CPU=24



USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
## this is only necessary if you want to add TEs to LTR families already named
##   for example, if you have another genome in this species that you want consistency for
##     it is not intentioned to compare across species, and also expects that these are all sequences of just the LTR, not the internal domains
LTRFILE=/group/jrigrp8/mstitzer/multiple_maize_assemblies/nested_ltr/families/B73V4.both_pseudo_AND_unplaced.allnested.5ltr.newNamewFam.fa


MEMLIM=96GB ## for step 2, nestedloop

MEMLIM=64GB ## for step 3, use entire med node
CPU=24









