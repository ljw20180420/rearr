#!/usr/bin/env python
import sys, fuzzysearch, numpy, os, Bio.Align
# _, csvfile, primer = sys.argv
csvfile, primer = "barcode/final_hgsgrna_libb_all_0811_NAA.csv", "TCAAGACCTAGCTAGCGAATT"
with os.popen(f'''tail -n+2 {csvfile} | cut -d',' -f12 | tr ACGT TGCA | rev ''', "r") as bd: # read barcodes
    barcodes = [barcode.rstrip() for barcode in bd]

min_seq_size = 20
min_score = -2
aligner = Bio.Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 0
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
aligner.target_end_gap_score = 0

# with os.popen("sed -n '2~4p' barcode/test3/wt2-g3n-1.R2.fq | head -n 400 | sort | uniq -c", "r") as pd:
#     for line in pd:
for line in sys.stdin:
    count, seq = line.strip().split()
    primer_alignment = aligner.align(primer, seq)[0]
    if primer_alignment.score >= min_score:
        primer_start, primer_end = primer_alignment.aligned[1][0][0], primer_alignment.aligned[1][-1][-1]
    else:
        _ = sys.stderr.write(f"{count}\t{seq}\tprimer_not_match\n")
        continue
    best_score = -numpy.inf
    for barcode in barcodes:
        barcode_alignment = aligner.align(barcode, seq[primer_end:])[0]
        if barcode_alignment.score > best_score:
            best_score = barcode_alignment.score
            best_barcode = barcode
            best_barcode_alignment = barcode_alignment
    if best_score < min_score:
        _ = sys.stderr.write(f"{count}\t{seq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\tbarcode_too_fuzzy\n")
        continue
    else:
        barcode_start, barcode_end = primer_end + best_barcode_alignment.aligned[1][0][0], primer_end + best_barcode_alignment.aligned[1][-1][-1]
        if barcode_end + 3 + min_seq_size > len(seq):
            sys.stderr.write(f"{count}\t{seq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{barcode_start}\t{barcode_end}\t{seq[barcode_start:barcode_end]}\tseq_downstream_to_barcode_too_short\n")
            continue
    _ = sys.stdout.write(f"{count}\t{seq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{barcode_start}\t{barcode_end}\t{seq[barcode_start:barcode_end]}\n")
