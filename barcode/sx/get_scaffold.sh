#!/bin/bash
# Usage: get_scaffold.sh <csvfile
# csvfile = 20bp + sgRNA(20bp) + scaffold(93bp) + target(44bp) + 3bp + RCbarcode(18bp) + RCprimer(21bp)
cut -d, -f2 | sed -r 's/^[ACGTN]+//; s/[ACGTN]+\s+$//' | dd conv=ucase 2>/dev/null