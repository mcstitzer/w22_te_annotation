bedtools getfasta -fi $GENOMENAME -bed ${GENOMENAME}.fa.up.gff3 -tab -fo ${GENOMENAME}.fa.up.tabout
bedtools getfasta -fi $GENOMENAME -bed ${GENOMENAME}.fa.down.gff3 -tab -fo ${GENOMENAME}.fa.down.tabout


