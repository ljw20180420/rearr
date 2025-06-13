#!/usr/bin/env python

import Bio.Align, sys, pysam, more_itertools, re


def coor2indel(refcoor, seqcoor, cut1, cut2):
    segs = []
    for i in range(len(refcoor) - 1):
        if refcoor[i + 1] > refcoor[i] and seqcoor[i + 1] > seqcoor[i]:
            segs.append([refcoor[i], refcoor[i + 1], seqcoor[i], seqcoor[i + 1]])
    indels = []
    for i in range(len(segs) - 1):
        indels.append(
            [segs[i][1], segs[i + 1][0] - cut1 + cut2, segs[i][3], segs[i + 1][2]]
        )
    if not indels:
        seqpos = segs[0][2] + cut1 - segs[0][0]
        indels.append([cut1, cut2, seqpos, seqpos])
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
    cut1 = int(sys.argv[2])
    cut2 = int(sys.argv[3])

    if program.lower() == "crispresso":
        sys.stdin.readline()
        for line in sys.stdin:
            align_seq, align_ref, _ = line.split("\t", 2)
            refcoor, seqcoor = Bio.Align.Alignment.parse_printed_alignment(
                [align_ref.encode(), align_seq.encode()]
            )[1]
            indels = coor2indel(refcoor, seqcoor, cut1, cut2)
            sys.stdout.write(f"{align_seq.replace('-', '')}")
            for indel in indels:
                sys.stdout.write(f"\t{indel[0]}\t{indel[1]}\t{indel[2]}\t{indel[3]}")
            sys.stdout.write("\n")

    if program.lower() == "crisprvariants":
        gquery, fc = "", ""
        for line in sys.stdin:
            query, flag, _, pos, _, cigar, _ = line.split("\t", 6)
            if (int(flag) / 4) % 2 or (int(flag) / 16) % 2:
                continue
            if query != gquery:
                sys.stdout.write(f"{fc}{query}")
                gquery, fc = query, "\n"
            refcoor, seqcoor = cigar2coor(cigar)
            indels = coor2indel(refcoor, seqcoor, cut1 - int(pos), cut2)
            for indel in indels:
                sys.stdout.write(
                    f"\t{indel[0] + int(pos)}\t{indel[1]}\t{indel[2]}\t{indel[3]}"
                )

    if program.lower() == "amplican":
        amplicon_start = int(sys.argv[4])
        primer_len = int(sys.argv[5])
        for _, align_ref, align_seq, _ in more_itertools.batched(sys.stdin, 4):
            refcoor, seqcoor = Bio.Align.Alignment.parse_printed_alignment(
                [align_ref.encode(), align_seq.encode()]
            )[1]
            indels = coor2indel(refcoor, seqcoor, cut1 - amplicon_start, cut2)
            for i in range(len(indels)):
                lc = "\t" if i < len(indels) - 1 else "\n"
                indel = indels[i]
                sys.stdout.write(
                    f"{indel[0] + amplicon_start}\t{indel[1]}\t{indel[2] - primer_len}\t{indel[3] - primer_len}{lc}"
                )

    if program.lower() == "amplicondivider":
        for line, _ in more_itertools.batched(sys.stdin, 2):
            query, _, ref, pos, _, cigar, _ = line.split("\t", 6)
            if cigar == "*":
                continue
            refcoor, seqcoor = cigar2coor(cigar)
            indels = coor2indel(refcoor, seqcoor, cut1 - int(pos), cut2)
            sys.stdout.write(f"{query}")
            for indel in indels:
                sys.stdout.write(
                    f"\t{indel[0] + int(pos)}\t{indel[1]}\t{indel[2]}\t{indel[3]}"
                )
            sys.stdout.write("\n")

    if program.lower() == "crispr-grant":
        gquery, fc = "", ""
        for line in sys.stdin:
            query, flag, _, pos, _, cigar, _ = line.split("\t", 6)
            if (int(flag) / 4) % 2 or (int(flag) / 16) % 2:
                continue
            if query != gquery:
                sys.stdout.write(f"{fc}{query}")
                gquery, fc = query, "\n"
            refcoor, seqcoor = cigar2coor(cigar)
            indels = coor2indel(refcoor, seqcoor, cut1 - int(pos), cut2)
            for indel in indels:
                sys.stdout.write(
                    f"\t{indel[0] + int(pos)}\t{indel[1]}\t{indel[2]}\t{indel[3]}"
                )

    if program.lower() == "zhangfeng":
        refstart = int(sys.argv[4])
        for line in sys.stdin:
            fields = line.rstrip().split()
            query = fields[0]
            ref_seq_coor = [int(s) for s in fields[1:]]
            indels = coor2indel(
                ref_seq_coor[: len(ref_seq_coor) // 2],
                ref_seq_coor[len(ref_seq_coor) // 2 :],
                cut1 - refstart,
                cut2,
            )
            sys.stdout.write(f"{query}")
            for indel in indels:
                sys.stdout.write(
                    f"\t{indel[0] + refstart}\t{indel[1]}\t{indel[2]}\t{indel[3]}"
                )
            sys.stdout.write("\n")

    if program.lower() == "selftarget":
        for line in sys.stdin:
            query, _, identifier, _ = line.split("\t", 3)
            try:
                (
                    _,
                    rpos1min,
                    rpos1max,
                    rpos2min,
                    rpos2max,
                    qpos1min,
                    qpos1max,
                    qpos2min,
                    qpos2max,
                ) = identifier.split("|")
            except Exception:
                continue
            if identifier == "-":
                continue
            sys.stdout.write(
                f"{query}\t{rpos1min}\t{int(rpos2min) - cut1 + cut2}\t{qpos1min}\t{qpos2min}\t{rpos1max}\t{int(rpos2max) - cut1 + cut2}\t{qpos1max}\t{qpos2max}\n"
            )
