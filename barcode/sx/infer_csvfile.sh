#!/bin/bash
# Usage: infer_csvfile.sh fastqR1
project_path="$(dirname $(realpath $0))/../.."
chip=$(gawk -F "-" '{print $(NF - 1)}' <<<$1 | head -c2)
ls ${project_path}/barcode/sx/csvfiles/*.csv | grep "${chip^^}"