#!/bin/bash
fastqfile=$1
csvfile=$2
primer=${3:-"TCAAGACCTAGCTAGCGAATT"}
sed -n '2~4p' $fastqfile | sort | uniq -c | python barcode/barcode.py $csvfile $primer 2> $fastqfile.$(basename $csvfile).not_find | sort -t $'\t' -k6,6 -k1,1nr > $fastqfile.$(basename $csvfile).barcode