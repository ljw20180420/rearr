#!/usr/bin/env python

import sys, numpy

def lrs_dis(lrs1, lrs2):
    s1s = [lrs1[0][0] + lrs1[1][0], lrs1[0][1] + lrs1[1][1]]
    s2s = [lrs2[0][0] + lrs2[1][0], lrs2[0][1] + lrs2[1][1]]
    if s1s[1] < s2s[0]:
        return [lrs1[0][1] - lrs2[0][0], lrs1[1][1] - lrs2[1][0]]
    elif s1s[0] > s2s[1]:
        return [lrs1[0][0] - lrs2[0][1], lrs1[1][0] - lrs2[1][1]]
    else:
        ld = (lrs1[0][0] - lrs2[0][0]) / 2 + (lrs2[1][0] - lrs1[1][0]) / 2
        return [ld, -ld]

def indel2lrs(indel, ref):
    lrs = [[indel[0]] * 2, [indel[1]] * 2]
    while lrs[0][0] - 1 > 0 and ref[lrs[0][0] - 1] == ref[lrs[1][0] - 1]:
        lrs[0][0], lrs[1][0] = lrs[0][0] - 1, lrs[1][0] - 1
    while lrs[1][1] < len(ref) and ref[lrs[0][1]] == ref[lrs[1][1]]:
        lrs[0][1], lrs[1][1] = lrs[0][1] + 1, lrs[1][1] + 1
    return lrs


def indel_dis(indel1, indel2, ref):
    lrs1, lrs2 = indel2lrs(indel1, ref), indel2lrs(indel2, ref)
    return lrs_dis(lrs1, lrs2) + [indel1[2] - indel2[2]]

ref = sys.argv[1]

for line in sys.stdin:
    mindis = numpy.inf
    query, l1, r1, m1, c2 = line.split("\t", 4)
    indel1 = [int(l1), int(r1), int(m1)]
    c2 = c2.split("\t")
    for i in range(0, len(c2), 3):
        indel2 = [int(c2[i]), int(c2[i + 1]), int(c2[i + 2])]
        disv = indel_dis(indel1, indel2, ref)
        dis = numpy.sum(numpy.abs(disv))
        if dis < mindis:
            mindis = dis
            mindisv = disv
            # minindel2 = indel2 # debug
    sys.stdout.write(f"{query}\t{mindisv[0]}\t{mindisv[1]}\t{mindisv[2]}\n")
    # sys.stdout.write(f"{query}\t{mindisv[0]}\t{mindisv[1]}\t{mindisv[2]}\t{indel1[0]}\t{indel1[1]}\t{indel1[2]}\t{minindel2[0]}\t{minindel2[1]}\t{minindel2[2]}\n") # debug