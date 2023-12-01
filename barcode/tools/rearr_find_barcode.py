#!/usr/bin/env python
import sys, pysam

for line in sys.stdin:
    query, _, barcode, _, _, CIGAR, _ = line.split("\t", 6)
    pAS = pysam.AlignedSegment()
    pAS.cigarstring = CIGAR
    sys.stdout.write(f"{int(query)+1}\t{barcode}\t{pAS.query_alignment_end}\n")
