#!/bin/bash

find_barcode()
{
    primer="TCAAGACCTAGCTAGCGAATT"
    minscore=30
    fastqfile=$1
    csvfile=$2

    cut -f1 "$fastqfile.count" | bowtie2 --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,$minscore -r -x "$csvfile.primer+barcode" -U - 2>/dev/null | samtools view | awk -F "\t" -v OFS="\t" '($2/4)%2!=0 || ($2/16)%2!=0{print > "'"$fastqfile.$(basename $csvfile).not_find"'"} ($2/4)%2==0 && ($2/16)%2==0{print}' | rearr_find_barcode.py | sort -k1,1 | join -t $'\t' -1 1 -2 1 <(cut -f1 "$fastqfile.count" | sed '=' | sed '1~2s/^/>s/' | cutadapt -a GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC - 2> /dev/null | sed '1~2d' | paste <(nl -w1 "$fastqfile.count") - | sort -k1,1) - | cut -f2- | sort -k4,4 | join -t $'\t' -1 4 -2 1 - <(sed -r 's/^>//; N; s/\n'"$primer"'/\t/' "$csvfile.primer+barcode.fa" | sort -k1,1) | cut -f2- | awk -v minscore=$minscore '$4+minscore>length($3){print >"'"$fastqfile.$(basename $csvfile).too_short"'"} $4+minscore<=length($3){print}' | sort -k5,5 -k2,2nr >"$fastqfile.$(basename $csvfile).barcode"
}

bowtie2genome=$1
getfastagenome=$2

while read fastqfile csvfile
do
    ( find_barcode "$fastqfile" "$csvfile"; rearr_barcode_align.py "$fastqfile.$(basename $csvfile).barcode" "$csvfile" "$bowtie2genome" "$getfastagenome" ) &
done

jobs
wait