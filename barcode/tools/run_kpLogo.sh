#!/bin/bash

# usage
# run_kpLogo.sh path_to_fastq method(thres|weight|back) threshold

fq1=$1
csvfile=$(infer_csvfile.sh "$fq1")
method=$2
threshold=${3:-40}
project_path="$(dirname $(realpath $0))/../.."

rearr_barcode_post_process.sh <"$fq1.table" | tail -n+2 | join -t $'\t' -1 1 -2 1 - <(rev "$csvfile" | cut -c23-40 | tr 'ACGT' 'TGCA' | paste - <(sed -r 's/[acgt]+.*$//' "$csvfile" | rev | cut -c1-20 | rev) | sort) | awk -F "\t" -v OFS="\t" -v fq1="$fq1" -v threshold="$threshold" '
    {
        if ($25 != sgRNA && sgRNA)
        {
            for (i = 0; i < inscount; ++i)
                print sgRNA
            for (i = 0; i < count; ++i)
                print sgRNA > fq1 ".kpLogo.total"
            print sgRNA, inscount / count, inscount, count > fq1 ".kpLogo.weight"
            if (count > threshold)
                print sgRNA, inscount / count, inscount, count > fq1 ".kpLogo.weight.thres"
            count = 0
            inscount = 0
        }
        sgRNA = $25
        count += $3
        if ($24 == "ins" || $24 == "indel")
            inscount += $3
    }
    END{
        for (i = 0; i < inscount; ++i)
            print sgRNA
        for (i = 0; i < count; ++i)
            print sgRNA > fq1 ".kpLogo.total"
        print sgRNA, inscount / count, inscount, count > fq1 ".kpLogo.weight"
        if (count > threshold)
            print sgRNA, inscount / count, inscount, count > fq1 ".kpLogo.weight.thres"
    }
' >"$fq1.kpLogo.insertion"

case $method in
    "thres")
        $project_path/barcode/kpLogo/bin/kpLogo "$fq1.kpLogo.weight.thres" -o "$fq1.kpLogo" -region 15,20 -weighted -k 1
    ;;
    "weighted"|"weight")
        $project_path/barcode/kpLogo/bin/kpLogo "$fq1.kpLogo.weight" -o "$fq1.kpLogo" -region 15,20 -weighted -k 1
    ;;
    "background"|"back"|"bg"|"bgfile")
        $project_path/barcode/kpLogo/bin/kpLogo "$fq1.kpLogo.insertion" -o "$fq1.kpLogo" -region 15,20 -bgfile "$fq1.kpLogo.total" -k 1
    ;;
esac
