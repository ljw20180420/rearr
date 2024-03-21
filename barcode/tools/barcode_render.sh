#!/bin/bash
# Usage
# barcode_render.sh fqR1 ext1up ext2up ref1len

fqR1=$1
ext1up=$2
ext2up=$3
ref1len=$4
algfile="${fqR1}.alg"
tablefile="${fqR1}.table"
sed -nr '/^[1-9]/,+2p; /^$/d' "${algfile}" >"${algfile}.${ext1up}.${ext2up}.${ref1len}"
cut -f1 --complement ${tablefile} >"${tablefile}.${ext1up}.${ext2up}.${ref1len}"
rearr_render.sh "${algfile}.${ext1up}.${ext2up}.${ref1len}"
rm "${tablefile}.${ext1up}.${ext2up}.${ref1len}"
rm "${algfile}.${ext1up}.${ext2up}.${ref1len}"