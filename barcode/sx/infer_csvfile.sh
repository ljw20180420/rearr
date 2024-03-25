#!/bin/bash
# Usage: infer_csvfile.sh fastqR1
project_path="$(dirname $(realpath $0))/../.."
chip=$("${project_path}/gawk-5.3.0/gawk" -F "-" '{print $(NF - 1)}' <<<$1 | head -c2)
ls ${project_path}/barcode/sx/csvfiles/*.csv | grep "${chip^^}"