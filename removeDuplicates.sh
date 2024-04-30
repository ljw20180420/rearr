#!/bin/bash
# Paste input fastq files side by side. Then remove and count duplicate rows.
# Usage: removeDuplicates.sh fastq1 fastq2 fastq3 ... >rmDupFile

paste "$@" | sed -n '2~4p' | sort | uniq -c | gawk '
    {
        for (i = 2; i <= NF; ++i)
            printf("%s\t", $i)
        print $1
    }
'