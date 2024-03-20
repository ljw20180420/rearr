#!/bin/bash

csvpath="$(dirname $(realpath $0))/csvfiles"
chip=$(awk -F "-" '{print $(NF - 1)}' <<<$1 | head -c2)
case ${chip^^} in
    G?)
        echo "$csvpath/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_${chip^^}.csv"
        ;;
    A?)
        echo "$csvpath/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_${chip^^}.csv"
        ;;
    *)
        echo "cannot infer csv file from fastq file" >&2
        exit
        ;;
esac