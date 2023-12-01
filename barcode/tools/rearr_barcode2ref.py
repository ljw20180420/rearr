#!/usr/bin/env python
import sys
_, extup_left, extdown_left, extup_right, extdown_right = sys.argv
extup_left, extdown_left, extup_right, extdown_right = int(extup_left), int(extdown_left), int(extup_right), int(extdown_right)
for line in sys.stdin:
    qname, flag, chr, pos, _, CIGAR, _ = line.split("\t", 6)
    flag, pos, tglen = int(flag), int(pos), int(CIGAR[:-1])
    tgstrand = "-" if (flag // 16) % 2 else "+"
    cut = pos - 1 + 16 + 6 if tgstrand == "+" else pos - 1 + tglen - 16 - 6
    rfstrand = "+" if tgstrand == "-" else "-" # the strand of target is constract to that of barcode and reference
    if rfstrand == "-":
        eul, edl, eur, edr = extdown_left, extup_left, extdown_right, extup_right
    else:
        eul, edl, eur, edr = extup_left, extdown_left, extup_right, extdown_right
    sys.stdout.write(f"{chr}\t{cut - eul}\t{cut + edl}\t{qname}\t.\t{rfstrand}\n{chr}\t{cut - eur}\t{cut + edr}\t{qname}\t.\t{rfstrand}\n")