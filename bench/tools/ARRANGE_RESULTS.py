#!/usr/bin/env python

import Bio.Align, sys, pysam, more_itertools, re

def coor2indel(refcoor, seqcoor):
    indels = []
    for i in range(len(refcoor) -1):
        if refcoor[i + 1] > refcoor[i] and seqcoor[i + 1] > seqcoor[i]:
            if indels:
                indels[-1][1] = refcoor[i]
            indels.append([refcoor[i + 1], None, 0])
        else:
            if not indels:
                indels.append([0, None, 0])
            if refcoor[i + 1] == refcoor[i]:
                indels[-1][2] += seqcoor[i + 1] - seqcoor[i]
    if not indels[-1][1]:
        indels[-1][1] = refcoor[-1]
    if indels[-1][0] == indels[-1][1] and indels[-1][2] == 0:
        indels.pop()
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
        
program = sys.argv[1]

if program.lower() == "crispresso":
    sys.stdin.readline()
    for line in sys.stdin:
        align_seq, align_ref, _ = line.split("\t", 2)
        refcoor, seqcoor = Bio.Align.Alignment.infer_coordinates([align_ref, align_seq])
        indels = coor2indel(refcoor, seqcoor)
        sys.stdout.write(f"{align_seq.replace('-', '')}")
        for indel in indels:
            sys.stdout.write(f"\t{indel[0]}\t{indel[1]}\t{indel[2]}")
        sys.stdout.write("\n")

if program.lower() == "crisprvariants":
    sgstart = int(sys.argv[2])
    for line in sys.stdin:
        query, cigar = line.split("\t")
        refcoor, seqcoor = cigar2coor(cigar)
        indels = coor2indel(refcoor, seqcoor)
        sys.stdout.write(f"{query}")
        for indel in indels:
            sys.stdout.write(f"\t{indel[0] + sgstart}\t{indel[1] + sgstart}\t{indel[2]}")
        sys.stdout.write("\n")

if program.lower() == "amplican":
    amplicon_start = int(sys.argv[2])
    for header, align_ref, align_seq, _ in more_itertools.batched(sys.stdin, 4):
        id = int(re.search("read_id:\s+(\d+)", header).group(1))
        refcoor, seqcoor = Bio.Align.Alignment.infer_coordinates([align_ref, align_seq])
        indels = coor2indel(refcoor, seqcoor)
        sys.stdout.write(f"@seq{id}")
        for indel in indels:
            sys.stdout.write(f"\t{amplicon_start + indel[0]}\t{amplicon_start + indel[1]}\t{indel[2]}")
        sys.stdout.write("\n")

if program.lower() == "amplicondivider":
    for line, _ in more_itertools.batched(sys.stdin, 2):
        query, _, ref, pos, _, cigar, _ = line.split('\t', 6)
        if ref == "*":
            continue
        refcoor, seqcoor = cigar2coor(cigar)
        indels = coor2indel(refcoor, seqcoor)
        sys.stdout.write(f"{query}")
        for indel in indels:
            sys.stdout.write(f"\t{indel[0] + int(pos)}\t{indel[1] + int(pos)}\t{indel[2]}")
        sys.stdout.write("\n")