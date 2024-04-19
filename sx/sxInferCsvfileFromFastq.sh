#!/bin/bash
# Usage: infer_csvfile.sh fastqR1 csvPath
chip=$(gawk -F "-" '{print $(NF - 1)}' <<<$1 | head -c2)
ls $2/*.csv | grep "${chip^^}"