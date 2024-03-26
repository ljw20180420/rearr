#!/bin/bash
project_path="$(dirname $(realpath $0))/../.."
fqR1="${project_path}/barcode/test/A2-g1n-3.fq"
fqR2="${project_path}/barcode/test/A2-g1n-3.R2.fq"
csvfile="$(${project_path}/barcode/sx/infer_csvfile.sh ${fqR1})"
spliter1="${csvfile}.adapter+sgRNA+scaffold"
spliter2="${csvfile}.primer+barcode"
sgRNAfile="${csvfile}.sgRNA"
ref12="${csvfile}.ref12"
# read -ep "path to genome reference:" genomeref

# $project_path/barcode/sx/index_spliter.sh $genomeref $(ls ${project_path}/barcode/sx/csvfiles/*.csv)
$project_path/barcode/tools/demultiplex.sh "${fqR1}" "${fqR2}" "${spliter1}" "${spliter2}" "${sgRNAfile}" "${ref12}" 100 30 30
$project_path/gawk-5.3.0/gawk -f $project_path/barcode/tools/barcode_align.AWK <"${fqR1}.demultiplex" -- $project_path "${fqR1}" 50 10 $(cat "${fqR1}.total")
"${project_path}/barcode/tools/run_kpLogo.sh" "$project_path/barcode/test/A2-g1n-3.fq" weight
"${project_path}/barcode/tools/trouble_shooting.r" "${fqR1}" -9999
"${project_path}/barcode/tools/barcode_render.sh" "${project_path}/barcode/test/A2-g1n-3.fq" 50 10 60