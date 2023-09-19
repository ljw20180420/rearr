#!/bin/bash
# while read barcode; do
#     grep --color $barcode A2_TEST.fq
# done < outbarcoede.txt
python barcode.py TCAAGACCTAGCTAGCGAATT outbarcoede.txt A2_TEST.fq | sort -t $'\t' -k7,7 | awk -F "\t" -v OFS="\t" '{if ($5 != cbarcode) {cbarcode = $5; if (NR > 1) print "----"}; print}' | csplit - --prefix='output/barcode' /----/ '{*}'


