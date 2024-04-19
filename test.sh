#!/bin/bash

for target in $@
do
    case $target in
    getSxCsvFileRef)
        read -ep "path to genome reference:" genomeref
        for csvfile in $(ls sx/csvfiles/*.csv)
        do
            getSxCsvFileRef.sh ${csvfile} ${genomeref} >${csvfile}.ref
        done
        ;;
    sxIndexSpliter)
        read -ep "path to genome reference:" genomeref
        sxIndexSpliter.sh genomeref $(ls sx/csvfiles/*.csv)
        ;;
    demultiplex)
        demultiplex.sh test/A2-g1n-3.fq.count sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.primer+barcode sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.adapter+sgRNA+scaffold 30 100 30
        ;;
    correct|rearrangement|Rearrangement)
        rearrangement <test/A2-g1n-3.fq.count.demultiplex 3<sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -5 | gawk -f correct_micro_homology.awk -- sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref NGG NGG >test/A2-g1n-3.fq.count.demultiplex.alg
        ;;
    esac
done