#!/bin/bash

primer="TCAAGACCTAGCTAGCGAATT"
csvfile=$1

tail -n+2 $csvfile | cut -d',' -f12 | tr ACGT TGCA | rev | sort | sed '=' | sed '1~2s/^/>BC_/; 2~2s/^/'"$primer"'/' >"$csvfile.primer+barcode.fa"
bowtie2-build -q "$csvfile.primer+barcode.fa" "$csvfile.primer+barcode"
