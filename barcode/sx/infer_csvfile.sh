#!/bin/bash
# Usage: infer_csvfile.sh fastqR1
csvpath="$(dirname $(realpath $0))/csvfiles"
chip=$(awk -F "-" '{print $(NF - 1)}' <<<$1 | head -c2)
ls $csvpath/*.csv | grep "${chip^^}"