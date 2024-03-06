#!/bin/bash
project_path="$(dirname $(realpath $0))/../.."
read -ep "path to genome reference:" genomeref
$project_path/barcode/tools/rearr_barcode_index.sh $genomeref
printf "%s\n" "$project_path/barcode/test/A2-g1n-3.fq" | "$project_path/barcode/tools/rearr_barcode_run.sh"
$project_path/barcode/tools/run_kpLogo.sh "$project_path/barcode/test/A2-g1n-3.fq" weight
$project_path/barcode/tools/trouble_shooting.r "$project_path/barcode/test/A2-g1n-3.fq" -9999
$project_path/barcode/tools/rearr_barcode_render.sh "$project_path/barcode/test/A2-g1n-3.fq"