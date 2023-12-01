#!/usr/bin/env python
import sys, fuzzysearch, numpy, os
_, countfile, csvfile, primer = sys.argv
# csvfile, primer = "barcode/final_hgsgrna_libb_all_0811_NAA.csv", "TCAAGACCTAGCTAGCGAATT"
with os.popen(f'''tail -n+2 {csvfile} | cut -d',' -f12 | tr ACGT TGCA | rev ''', "r") as bd: # read barcodes
    barcodes = [barcode.rstrip() for barcode in bd]
    barcode_mat = numpy.array([list(barcode) for barcode in barcodes])
barcode_size = len(barcodes[0])
min_barcode_size = barcode_size - 2
shift_len = 2 # shift length to check barcode shift
max_l_dist = 2
min_seq_size = 20

with open(countfile, "r") as cd, os.popen(f''' cut -f1 {countfile} | sed '=' | sed '1~2s/^/>s/' | cutadapt -a GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC - 2> /dev/null | sed '1~2d' ''', "r") as ad:
    for line in cd:
        seq, count = line.rstrip().split()
        naseq = ad.readline().rstrip()
        primer_start = seq.find(primer)
        if primer_start != -1:
            primer_end = primer_start + len(primer)
        else:
            primer_matches = fuzzysearch.find_near_matches(primer, seq, max_l_dist = max_l_dist)
            if not primer_matches:
                _ = sys.stderr.write(f"{count}\t{seq}\tprimer_not_match\n")
                continue
            primer_match = min(primer_matches, key=lambda match: match.dist)
            primer_start, primer_end = primer_match.start, primer_match.end
        barcode = seq[primer_end:min(primer_end+barcode_size, len(seq))]
        if len(barcode) < min_barcode_size:
            sys.stderr.write(f"{count}\t{seq}\t{naseq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\tbarcode_too_short\n")
            continue
        if barcode in barcodes:
            barcode_start, barcode_end = primer_end, primer_end + barcode_size
        else:
            permeable_shift_len = min(shift_len, len(seq) - primer_end - len(barcode))
            barcode_shifts = numpy.array([list(seq[primer_end + shif:primer_end + len(barcode) + shif]) for shif in range(permeable_shift_len + 1)]).reshape(permeable_shift_len + 1, 1, -1)
            barcode = barcodes[numpy.argmin(numpy.min(numpy.sum(barcode_mat[:,:len(barcode)] != barcode_shifts, axis=2), axis=0))]
            barcode_matches = fuzzysearch.find_near_matches(barcode, seq[primer_end:], max_l_dist = max_l_dist)
            if not barcode_matches:
                _ = sys.stderr.write(f"{count}\t{seq}\t{naseq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\tbarcode_too_fuzzy\n")
                continue                        
            barcode_match = min(barcode_matches, key=lambda match: match.dist)
            barcode_start, barcode_end = primer_end + barcode_match.start, primer_end + barcode_match.end
        if barcode_end + 3 + min_seq_size > len(naseq):
            sys.stderr.write(f"{count}\t{seq}\t{naseq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{barcode_start}\t{barcode_end}\t{seq[barcode_start:barcode_end]}\tseq_downstream_to_barcode_too_short\n")
            continue
        _ = sys.stdout.write(f"{count}\t{seq}\t{naseq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{barcode_start}\t{barcode_end}\t{seq[barcode_start:barcode_end]}\n")
