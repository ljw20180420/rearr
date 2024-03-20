#!/bin/bash
# Usage: get_barcode.sh <csvfile
# csvfile = 20bp + sgRNA(20bp) + scaffold(83/93bp) + target(44bp) + 3bp + RCbarcode(18bp) + RCprimer(21bp)
rev | sed -r 's/^\s+//' | cut -c22-39 | dd conv=ucase 2>/dev/null | tr 'ACGT' 'TGCA'