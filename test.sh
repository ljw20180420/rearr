#!/bin/bash

for target in $@
do
    case $target in
    removeDuplicates)
        removeDuplicates.sh test/A2-g1n-3.R2.fq test/A2-g1n-3.fq >test/A2-g1n-3.fq.count
        ;;
    sxCutR2AdapterFilterCumulate)
        sxCutR2AdapterFilterCumulate.sh test/A2-g1n-3.fq.count.demultiplex 30 >test/A2-g1n-3.fq.count.demultiplex.post
        ;;
    getSxCsvFileRef)
        read -ep "path to genome reference:" genomeref
        for csvfile in $(ls sx/csvfiles/*.csv)
        do
            getSxCsvFileRef.sh ${csvfile} ${genomeref} >${csvfile}.ref
        done
        ;;
    sxExtractSpliter)
        read -ep "path to genome reference:" genomeref
        sxExtractSpliter.sh genomeref $(ls sx/csvfiles/*.csv)
        ;;
    demultiplex)
        demultiplex.sh test/A2-g1n-3.fq.count sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.primer+barcode sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.adapter+sgRNA+scaffold 30 100 >test/A2-g1n-3.fq.count.demultiplex
        ;;
    correct|rearrangement|Rearrangement)
        rearrangement <test/A2-g1n-3.fq.count.demultiplex.post 3<sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -5 | gawk -f correct_micro_homology.awk -- sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref NGG NGG >test/A2-g1n-3.fq.count.demultiplex.post.alg
        ;;
    esac
done