import sys, more_itertools, fuzzysearch, numpy

# primer = "TCAAGACCTAGCTAGCGAATT"
# barcodefile = "outbarcoede.txt"
# fastqfile = "A2_TEST.fq"
primer = sys.argv[1]
barcodefile = sys.argv[2]
fastqfile = sys.argv[3]
with open(barcodefile) as bd: # read barcodes
    barcodes = [barcode.rstrip() for barcode in bd]
    barcode_mat = numpy.array([list(barcode) for barcode in barcodes])
barcode_size = len(barcodes[0])
with open(fastqfile, "r") as fq:
    for query_name, seq, _, qual in more_itertools.batched(fq, 4):
        query_name, seq, qual = query_name.rstrip(), seq.rstrip(), qual.rstrip()
        primer_start = seq.find(primer)
        if primer_start != -1:
            primer_end = primer_start + len(primer)
        else:
            primer_matches = fuzzysearch.find_near_matches(primer, seq, max_l_dist=2)
            if not primer_matches:
                _ = sys.stdout.write(f"{query_name}\t{seq}\t{qual}\t*\n")
                continue
            primer_match = min(primer_matches, key=lambda match: match.dist)
            primer_start, primer_end = primer_match.start, primer_match.end
        barcode = seq[primer_end:primer_end+barcode_size]
        if barcode in barcodes:
            fuzzy_barcode_score, barcode_start, barcode_end = 100, primer_end, primer_end + barcode_size
        else:
            # continue
            barcode = barcodes[numpy.argmin(numpy.sum(barcode_mat == list(barcode), axis=1))]
            barcode_matches = fuzzysearch.find_near_matches(barcode, seq[primer_end:], max_l_dist=2)
            if not barcode_matches:
                _ = sys.stdout.write(f"{query_name}\t{seq}\t{qual}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{fuzzy_barcode_score}\t*\n")
                continue                           
            barcode_match = min(barcode_matches, key=lambda match: match.dist)
            barcode_start, barcode_end = primer_end + barcode_match.start, primer_end + barcode_match.end
        _ = sys.stdout.write(f"{query_name}\t{seq}\t{qual}\t{primer_start}\t{primer_end}\t{seq[primer_start:primer_end]}\t{barcode}\t{fuzzy_barcode_score}\t{barcode_start}\t{barcode_end}\t{seq[barcode_start:barcode_end]}\n")