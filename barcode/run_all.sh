#!/bin/bash

# while read fastqfile csvfile
# do
#     ( barcode/barcode.sh barcode/$fastqfile barcode/$csvfile; python barcode/align.py barcode/$fastqfile.$csvfile.barcode ) &
# done < <(echo "A2_TEST.fq" "final_hgsgrna_libb_all_0811-NGG.csv")
# wait

#this is the example to run many fq csv pairs
while read fastqfile csvfile
do
    ( barcode/barcode.sh barcode/$fastqfile barcode/$csvfile; python barcode/align.py barcode/$fastqfile.$csvfile.barcode ) &
done < <(echo "A2_TEST.fq" "final_hgsgrna_libb_all_0811-NGG.csv"$'\n' \
            "A2_TEST2.fq" "final_hgsgrna_libb_all_0811-NGG.csv")
wait