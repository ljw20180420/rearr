#!/usr/bin/env python
import sys, pysam

for line in sys.stdin:
    flag, barcode, CIGAR = line.split("\t")
    pAS = pysam.AlignedSegment()
    pAS.cigarstring = CIGAR
    sys.stdout.write(f"{flag}\t{barcode}\t{pAS.query_alignment_end}\n")
