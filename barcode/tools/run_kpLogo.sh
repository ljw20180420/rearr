#!/bin/bash

fq1=$1
csvfile=$(infer_csvfile.sh "$fq1")

rearr_barcode_post_process.sh <"$fq1.table" | tail -n+2 | join -t $'\t' -1 1 -2 1 - <(rev "$csvfile" | cut -c23-40 | tr 'ACGT' 'TGCA' | paste - <(sed -r 's/[acgt]+.*$//' "$csvfile" | rev | cut -c1-20 | rev) | sort) | awk -F "\t" -v OFS="\t" '
    $23 == "ins" || $23 == "indel"{
        if ($24 != sgRNA && sgRNA)
        {
            print sgRNA, count
            count = 0
        }
        sgRNA = $24
        count += $3
    }
    END{
        print sgRNA, count
    }
' >"$fq1.kpLogo"

barcode/kpLogo/bin/kpLogo "$fq1.kpLogo" -o "$fq1.kpLogo" -region 15,20