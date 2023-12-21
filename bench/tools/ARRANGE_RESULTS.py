#!/usr/bin/env python

import Bio.Align, sys, pysam, more_itertools, re

def coor2indel(refcoor, seqcoor, cut):
    segs = []
    for i in range(len(refcoor) - 1):
        if refcoor[i + 1] > refcoor[i] and seqcoor[i + 1] > seqcoor[i]:
            segs.append([refcoor[i], refcoor[i+1], seqcoor[i], seqcoor[i+1]])
    indels = []
    for i in range(len(segs) - 1):
        indels.append([segs[i][1], segs[i+1][0], segs[i][3], segs[i+1][2]])
    if not indels:
        seqpos = segs[0][2] + cut - segs[0][0]
        indels.append([cut, cut, seqpos, seqpos])
    return indels

def cigar2coor(cigar):
    tuplestring = "MIDNSHP=XB"
    refcoor, seqcoor = [0], [0]
    pAS = pysam.AlignedSegment()
    pAS.cigarstring = cigar
    for op, sz in pAS.cigartuples:
        refnc = refcoor[-1] + sz if tuplestring[op] in "M=XDN" else refcoor[-1]
        seqnc = seqcoor[-1] + sz if tuplestring[op] in "M=XISH" else seqcoor[-1]
        refcoor.append(refnc)
        seqcoor.append(seqnc)
    return refcoor, seqcoor

if __name__ == "__main__":
    program = sys.argv[1]
    cut = int(sys.argv[2])

    if program.lower() == "crispresso":
        sys.stdin.readline()
        for line in sys.stdin:
            align_seq, align_ref, _ = line.split("\t", 2)
            refcoor, seqcoor = Bio.Align.Alignment.infer_coordinates([align_ref, align_seq])
            indels = coor2indel(refcoor, seqcoor, cut)
            sys.stdout.write(f"{align_seq.replace('-', '')}")
            for indel in indels:
                sys.stdout.write(f"\t{indel[0]}\t{indel[1]}\t{indel[2]}\t{indel[3]}")
            sys.stdout.write("\n")

    if program.lower() == "crisprvariants":
        gquery = ""
        for i, line in enumerate(sys.stdin):
            query, flag, _, pos, _, cigar, _ = line.split('\t', 6)
            if (int(flag) / 4) % 2 or (int(flag) / 16) % 2:
                continue
            if query != gquery:
                fc = "" if i == 0 else "\n"
                sys.stdout.write(f"{fc}{query}")
                gquery = query
            refcoor, seqcoor = cigar2coor(cigar)
            indels = coor2indel(refcoor, seqcoor, cut - int(pos))
            for indel in indels:
                sys.stdout.write(f"\t{indel[0] + int(pos)}\t{indel[1] + int(pos)}\t{indel[2]}\t{indel[3]}")

    if program.lower() == "amplican":
        amplicon_start = int(sys.argv[3])
        primer_len = int(sys.argv[4])
        for _, align_ref, align_seq, _ in more_itertools.batched(sys.stdin, 4):
            refcoor, seqcoor = Bio.Align.Alignment.infer_coordinates([align_ref, align_seq])
            indels = coor2indel(refcoor, seqcoor, cut - amplicon_start)
            for i in range(len(indels)):
                lc = "\t" if i < len(indels) - 1 else "\n"
                indel = indels[i]
                sys.stdout.write(f"{amplicon_start + indel[0]}\t{amplicon_start + indel[1]}\t{indel[2] - primer_len}\t{indel[3] - primer_len}{lc}")

    if program.lower() == "amplicondivider":
        for line, _ in more_itertools.batched(sys.stdin, 2):
            query, _, ref, pos, _, cigar, _ = line.split('\t', 6)
            if cigar == "*":
                continue
            refcoor, seqcoor = cigar2coor(cigar)
            indels = coor2indel(refcoor, seqcoor, cut - int(pos))
            sys.stdout.write(f"{query}")
            for indel in indels:
                sys.stdout.write(f"\t{indel[0] + int(pos)}\t{indel[1] + int(pos)}\t{indel[2]}\t{indel[3]}")
            sys.stdout.write("\n")

    if program.lower() == "crispr-grant":
        for line in sys.stdin:
            query, flag, _, pos, _, cigar, _ = line.split('\t', 6)
            if (int(flag) / 4) % 2 or (int(flag) / 16) % 2:
                continue
            refcoor, seqcoor = cigar2coor(cigar)
            indels = coor2indel(refcoor, seqcoor, cut - int(pos))
            sys.stdout.write(f"{query}")
            for indel in indels:
                sys.stdout.write(f"\t{indel[0] + int(pos)}\t{indel[1] + int(pos)}\t{indel[2]}\t{indel[3]}")
            sys.stdout.write("\n")

    if program.lower() == "zhangfeng":
        refstart = int(sys.argv[3])
        for line in sys.stdin:
            fields = line.rstrip().split()
            query = fields[0]
            ref_seq_coor = [int(s) for s in fields[1:]]
            # sys.stderr.write(f"{ref_seq_coor}\n")
            # sys.stderr.write(f"{len(ref_seq_coor)/2}\n")
            indels = coor2indel(ref_seq_coor[:len(ref_seq_coor) // 2], ref_seq_coor[len(ref_seq_coor) // 2:], cut - refstart )
            sys.stdout.write(f"{query}")
            for indel in indels:
                sys.stdout.write(f"\t{indel[0] + refstart}\t{indel[1] + refstart}\t{indel[2]}\t{indel[3]}")
            sys.stdout.write("\n")

    if program.lower() == "selftarget":
        for line in sys.stdin:
            query, _, identifier_uw2_dw1, _ = line.split("\t", 3)
            try:
                identifier, uw2m1, dw1 = identifier_uw2_dw1.split("|")
            except Exception:
                continue
            if identifier == "-":
                continue
            lpos = int(re.search("L(-?\d+)", identifier).group(1)) + 1
            rpos = int(re.search("R(-?\d+)", identifier).group(1))
            sys.stdout.write(f"{query}\t{lpos + cut}\t{rpos + cut}\t{int(uw2m1) + 1}\t{int(dw1)}\n")