#!/bin/bash

fq1=$1
csvfile=$(infer_csvfile.sh "$fq1")

rearr_barcode_post_process.sh <"$fq1.table" | tail -n+2 | join -t $'\t' -1 1 -2 1 - <(rev "$csvfile" | cut -c23-40 | tr 'ACGT' 'TGCA' | paste - <(sed -r 's/[acgt]+.*$//' "$csvfile" | rev | cut -c1-20 | rev) | sort) | awk -F "\t" -v OFS="\t" -v fq1="$fq1" '
    {
        if ($24 != sgRNA && sgRNA)
        {
            for (i = 0; i < inscount; ++i)
                print sgRNA
            for (i = 0; i < count; ++i)
                print sgRNA > fq1 ".kpLogo.total"
            count = 0
            inscount = 0
        }
        sgRNA = $24
        count += $3
        if ($23 == "ins" || $23 == "indel")
            inscount += $3
    }
    END{
        for (i = 0; i < inscount; ++i)
            print sgRNA
        for (i = 0; i < count; ++i)
            print sgRNA > fq1 ".kpLogo.total"
    }
' >"$fq1.kpLogo.insertion"

barcode/kpLogo/bin/kpLogo "$fq1.kpLogo.insertion" -o "$fq1.kpLogo" -region 15,20 -bgfile "$fq1.kpLogo.total" -k 1

# cut -c15 "$fq1.kpLogo.insertion" | sort | uniq -c