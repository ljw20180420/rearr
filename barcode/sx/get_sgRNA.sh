#!/bin/bash
# Usage: get_sgRNA.sh <csvfile
# csvfile = 20bp + sgRNA(20bp) + scaffold(83/93bp) + target(44bp) + 3bp + RCbarcode(18bp) + RCprimer(21bp)
cut -d, -f2 | cut -c21-40 | dd conv=ucase 2>/dev/null