#!/bin/bash

chip=$(cut -d'-' -f2 <<<$1 | head -c2)
case ${chip^^} in
    G?)
        echo "barcode/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_${chip^^}.csv"
        ;;
    A?)
        echo "barcode/csvfiles/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_${chip^^}.csv"
        ;;
    *)
        echo "cannot infer csv file from fastq file" >&2
        exit
        ;;
esac