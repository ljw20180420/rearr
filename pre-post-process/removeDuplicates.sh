#!/bin/bash
# Paste input fastq files side by side. Then remove and count duplicate rows.
# Usage: removeDuplicates.sh fastq1(.gz|.zip) fastq2(.gz|.zip) fastq3(.gz|.zip) ... >rmDupFile

fqlist=""
for fq in "$@"
do
    if (file $fq | grep -q compressed)
    then
        fqlist="$fqlist <(zcat $fq)"
    else
        fqlist="$fqlist $fq"
    fi
done

eval paste $fqlist | sed -n '2~4p' | sort | uniq -c | gawk '
    {
        for (i = 2; i <= NF; ++i)
            printf("%s\t", $i)
        print $1
    }
'