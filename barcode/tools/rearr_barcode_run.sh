#!/bin/bash

find_barcode()
{
    sed -n '2~4p' "$1" | sort | uniq -c | awk -v OFS="\t" '{print $2, $1}' >"$1.count"
    rearr_find_barcode.py "$1.count" "$2" "TCAAGACCTAGCTAGCGAATT" 2> "$1.$(basename $2).not_find" | sort -t $'\t' -k7,7 -k1,1nr >"$1.$(basename $2).barcode"
}

bowtie2genome=$1
getfastagenome=$2

while read fastqfile csvfile
do
    ( find_barcode "$fastqfile" "$csvfile"; rearr_barcode_align.py "$fastqfile.$(basename $csvfile).barcode" "$csvfile" "$bowtie2genome" "$getfastagenome" ) &
done

jobs
wait