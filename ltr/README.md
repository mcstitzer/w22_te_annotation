## this replicates the nested_ltr approach 

### details at [maize_v4_te_annotation](https://github.com/mcstitzer/maize_v4_TE_annotation/tree/master/ltr)


## To predict structural LTRs:

- download tRNA and GyDb HMMs using ```get_tRNA_hmm_dbs.sh```, which are needed for ltrdigest

- ```run_ltrharvest.sh``` runs [ltrharvest](http://www.zbh.uni-hamburg.de/fileadmin/gi/LTRharvest/ltrharvestman.pdf) and [ltrdigest](http://www.zbh.uni-hamburg.de/fileadmin/gi/LTRdigest/ltrdigestman.pdf) on the genome

- but LTR TEs are nested, so we need to remove these copies and rerun. This is done in ```mask_subtract```


## To cluster into families

We use the 808080 rule ([Wicker et al., 2007](http://www.nature.com/nrg/journal/v8/n12/full/nrg2165.html)) to cluster. 

In the ```families``` directory:

- `cluster_families_LTR.sh` runs the initial allvsall vsearch of 5' LTRs and then goes through silix clustering of sequences

- `cluster_families_MTEC.sh` finds the closest match of each complete element to a [MTEC](http://maizetedb.org) consensus

## Combining all outputs

We assign superfamily and family names to each copy, return all sequences to the coordinates of the original genome, and output a gff3 of all LTR copies.

- `Rscript combine_all_gffs_and_assign_fams.R` 

	

