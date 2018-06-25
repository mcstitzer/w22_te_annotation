## Scripts used in TE annotation of *Zea mays* inbred line W22 for [Springer et al., 2018]()

The main differences from the TE annotation of the inbred line B73 in [Jiao et al., 2017](https://www.nature.com/articles/nature22971), [maize_v4_te_annotation](https://github.com/mcstitzer/maize_v4_TE_annotation) are:

- We add individual copies to existing families using the [80-80-80 rule](https://www.nature.com/articles/nrg2165).
	- in the interest of file size, fastas of existing families are not included in this repository. They can be generated from the B73v4 TE gff [here](ftp://ftp.gramene.org/pub/gramene/CURRENT_RELEASE/gff3/zea_mays/repeat_annotation/) using bedtools.
	- TEs that lack a family are clustered into new (W22-specific) families.

- We do not run detectMITE on the W22 assembly. Instead, we use detectMITE results from B73 as queries for target_mTEA searches for TIRs and TSDs.


Each order of TE was run separately. 
To make an overall gff, where TEs from different orders that overlap are filtered out, we first run `combine_all_LTRs.R`, then `combine_all_TEs.R` to generate the final gff3  (`W22.allTE.gff3`).	

A caveat: these scripts are designed for the FARM cluster at UC Davis, running SLURM. Additionally, scripts have many hard coded paths and dependencies. 
