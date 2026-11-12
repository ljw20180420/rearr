#!/usr/bin/env python

import sys

import Bio.Align

if __name__ == "__main__":
    program = sys.argv[1]
    cut1 = int(sys.argv[2])
    cut2 = int(sys.argv[3])

    if program.lower() == "crispr-a":
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

    if program.lower() == "siq":
        for line in sys.stdin:
            seq_name, rpos1, rpos2, qpos1, qpos2 = line.split("\t")
            rpos1, rpos2, qpos1, qpos2 = int(rpos1), int(rpos2), int(qpos1), int(qpos2)
            sys.stdout.write(
                f"{seq_name}\t{rpos1}\t{rpos2 - cut1 + cut2}\t{qpos1}\t{qpos2}\n"
            )
