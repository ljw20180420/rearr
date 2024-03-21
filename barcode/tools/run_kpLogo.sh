#!/bin/bash

# usage
# run_kpLogo.sh path_to_fastq spliter2file sgRNAfile method(thres|weight|back) threshold

project_path="$(dirname $(realpath $0))/../.."
fqR1=$1
spliter2file=$2
sgRNAfile=$3
method=$4
threshold=${5:-40}

"${project_path}/barcode/tools/barcode_post_process.sh" <"${fqR1}.table" | tail -n+2 | join -t $'\t' -1 1 -2 1 - <(paste ${spliter2file} ${sgRNAfile} | sort) | awk -F "\t" -v OFS="\t" -v fqR1="${fqR1}" -v threshold="$threshold" '
{
    if ($25 != sgRNA && sgRNA)
    {
        for (i = 0; i < inscount; ++i)
            print sgRNA
        for (i = 0; i < count; ++i)
            print sgRNA > fqR1 ".kpLogo.total"
        print sgRNA, inscount / count, inscount + 0, count > fqR1 ".kpLogo.weight"
        if (count > threshold)
            print sgRNA, inscount / count, inscount + 0, count > fqR1 ".kpLogo.weight.thres"
        count = 0
        inscount = 0
    }
    sgRNA = $25
    count += $3
    if ($24 == "ins" || $24== "indel")
        inscount += $3
}
END{
    for (i = 0; i < inscount; ++i)
        print sgRNA
    for (i = 0; i < count; ++i)
        print sgRNA > fqR1 ".kpLogo.total"
    print sgRNA, inscount / count, inscount + 0, count > fqR1 ".kpLogo.weight"
    if (count > threshold)
        print sgRNA, inscount / count, inscount + 0, count > fqR1 ".kpLogo.weight.thres"
}
' >"${fqR1}.kpLogo.insertion"

case $method in
    "thres")
        $project_path/barcode/kpLogo-1.1/bin/kpLogo "${fqR1}.kpLogo.weight.thres" -o "${fqR1}.kpLogo" -region 15,20 -weighted -k 1
    ;;
    "weighted"|"weight")
        $project_path/barcode/kpLogo-1.1/bin/kpLogo "${fqR1}.kpLogo.weight" -o "${fqR1}.kpLogo" -region 15,20 -weighted -k 1
    ;;
    "background"|"back"|"bg"|"bgfile")
        $project_path/barcode/kpLogo-1.1/bin/kpLogo "${fqR1}.kpLogo.insertion" -o "${fqR1}.kpLogo" -region 15,20 -bgfile "${fqR1}.kpLogo.total" -k 1
    ;;
esac
