These could likely be condensed to one script.

first, ```run_repeatmasker.sh``` uses the LTR annotation to identify homologous LTR regions in the genome assembly.

Then, ```repeatmasker_gff_to_flanking_gff.R``` generates and writes to file two files - a gff3 for the upstream TSD and another for the downstream TSD

Next, ```get_tsd_tabout.sh``` gets the sequence of the TSD using bedtools and writes to tabout

Then, ```flanking_gff_to_SOLO.R'``` finds those repeatmasker hits that have identical TSDs and filters the overall gff if they do not span existing LTR boundaries.
