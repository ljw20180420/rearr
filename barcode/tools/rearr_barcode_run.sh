#!/bin/bash

find_barcode()
{
    primer="TCAAGACCTAGCTAGCGAATT"
    minscore=30
    fastqfile=$1
    csvfile=$2

    sed -n '2~4p' "$fastqfile" | sort | uniq -c | awk -v OFS="\t" '{print $2, $1}' >"$fastqfile.count"
    cut -f1 "$fastqfile.count" | bowtie2 --quiet --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,$minscore -r -x "$csvfile.primer+barcode" -U - 2>/dev/null | samtools view | cut -f2,3,6 | $(which python) -c '
import sys, pysam
for line in sys.stdin:
    flag, barcode, CIGAR = line.split("\t")
    pAS = pysam.AlignedSegment()
    pAS.cigarstring = CIGAR
    sys.stdout.write(f"{flag}\t{barcode}\t{pAS.query_alignment_end}\n")
    ' | awk -F "\t" -v csvfile=$csvfile -v primer=$primer -v OFS="\t" 'BEGIN{pbfa = csvfile ".primer+barcode.fa"; while (getline bn <pbfa) {getline bc <pbfa; sub("^" primer, "", bc); bn2bc[substr(bn, 2)] = bc}} {print $1, bn2bc[$2], $3}' | paste "$fastqfile.count" <(cut -f1 "$fastqfile.count" | sed '=' | sed '1~2s/^/>s/' | cutadapt -a GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC - 2> /dev/null | sed '1~2d') - | awk -F "\t" -v OFS="\t" -v minscore=$minscore '{if (($4/4)%2==1 || ($4/16)%2==1) print $1, $2 > "'"$fastqfile.$(basename $csvfile).not_find"'"; else if ($6+minscore>length($3)) print $1, $2 > "'"$fastqfile.$(basename $csvfile).too_short"'"; else print $1, $2, $3, $6, $5}' | sort -k5,5 -k2,2nr >"$fastqfile.$(basename $csvfile).barcode"
}

bowtie2genome=$1
getfastagenome=$2

while read fastqfile csvfile
do
    ( find_barcode "$fastqfile" "$csvfile"; rearr_barcode_align.py "$fastqfile.$(basename $csvfile).barcode" "$csvfile" "$bowtie2genome" "$getfastagenome" ) &
done

jobs
wait


