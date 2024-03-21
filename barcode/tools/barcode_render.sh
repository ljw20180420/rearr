#!/bin/bash
# Usage
# barcode_render.sh fqfile

algfile=$1.alg
tablefile=$1.table
sed -nr '/^[1-9]/,+2p; /^$/d' $algfile >$algfile.50.10.60
cut -f1 --complement $tablefile >$tablefile.50.10.60
rearr_render.sh $algfile.50.10.60
rm $tablefile.50.10.60
rm $algfile.50.10.60