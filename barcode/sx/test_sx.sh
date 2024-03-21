#!/bin/bash
project_path="$(dirname $(realpath $0))/../.."
fqR1="${project_path}/barcode/test/A2-g1n-3.fq"
fqR2="${project_path}/barcode/test/A2-g1n-3.R2.fq"
csvfile="$($project_path/barcode/sx/infer_csvfile.sh $fqR1)"
spliter1="$csvfile.adapter+sgRNA+scaffold"
spliter2="$csvfile.primer+barcode"
ref12="$csvfile.ref12"
# read -ep "path to genome reference:" genomeref

# $project_path/barcode/sx/index_spliter.sh $genomeref
# $project_path/barcode/tools/demultiplex.sh "${fqR1}" "${fqR2}" "${spliter1}" "${spliter2}" 100 30 $ref12 30 >"${fqR1}.barcode"
# $project_path/barcode/tools/barcode_align.sh "${fqR1}" 50 10 $(cat "${fqR1}.total") <"${fqR1}.barcode" >"${fqR1}.alg" 3>"${fqR1}.table"
"${project_path}/barcode/tools/run_kpLogo.sh" "$project_path/barcode/test/A2-g1n-3.fq" <("${project_path}/barcode/sx/get_barcode.sh" <${csvfile}) <("${project_path}/barcode/sx/get_sgRNA.sh" <${csvfile}) weight
# $project_path/barcode/tools/trouble_shooting.r "$project_path/barcode/test/A2-g1n-3.fq" -9999
# $project_path/barcode/tools/rearr_render.sh "$project_path/barcode/test/A2-g1n-3.fq"