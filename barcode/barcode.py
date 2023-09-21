#!/home/xiaoge/miniconda3/envs/crispr/bin/python
import sys, fuzzysearch, numpy, os
csvfile = sys.argv[1]
primer = sys.argv[2]
with os.popen(f'''tail -n+2 {csvfile} | cut -d',' -f12 | tr ACGT TGCA | rev ''', "r") as bd: # read barcodes
    barcodes = [barcode.rstrip() for barcode in bd]
    barcode_mat = numpy.array([list(barcode) for barcode in barcodes])
barcode_size = len(barcodes[0])
shift_len = 3 # shift length to check barcode shift

for line in sys.stdin:
    count, seq = line.strip().split()
    primer_start = seq.find(primer)
    if primer_start != -1:
        primer_end = primer_start + len(primer)
    else:
        primer_matches = fuzzysearch.find_near_matches(primer, seq, max_l_dist=2)
        if not primer_matches:
            _ = sys.stderr.write(f"{count}\t{seq}\t*\n")
            continue
        primer_match = min(primer_matches, key=lambda match: match.dist)
        primer_start, primer_end = primer_match.start, primer_match.end
    barcode = seq[primer_end:primer_end+barcode_size]
    if barcode in barcodes:
        barcode_start, barcode_end = primer_end, primer_end + barcode_size
    else:
        barcode_shifts = numpy.array([list(seq[primer_end + shif:primer_end+barcode_size + shif]) for shif in range(shift_len)]).reshape(shift_len, 1, -1)
        barcode = barcodes[numpy.argmin(numpy.min(numpy.sum(barcode_mat == barcode_shifts, axis=1), axis=1))]
        barcode_matches = fuzzysearch.find_near_matches(barcode, seq[primer_end:], max_l_dist=2)
        if not barcode_matches:
            _ = sys.stderr.write(f"{count}\t{seq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t*\n")
            continue                        
        barcode_match = min(barcode_matches, key=lambda match: match.dist)
        barcode_start, barcode_end = primer_end + barcode_match.start, primer_end + barcode_match.end
    _ = sys.stdout.write(f"{count}\t{seq}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{barcode_start}\t{barcode_end}\t{seq[barcode_start:barcode_end]}\n")
