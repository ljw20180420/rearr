#!/usr/bin/env python

import sys, numpy

def NDsegment_distance(indel1, indel2, range1, range2):
    range = [range1[0] - range2[1], range1[1] - range2[0]]
    target = numpy.median(indel2 - indel1)
    actual = range[1] if target > range[1] else range[0] if target < range[0] else target
    return indel1 + actual - indel2

def get_range(indel, ref):
    range = [0, 0]
    while indel[0] + range[0] - 1 >= 0 and indel[1] + range[0] - 1 >= 0 and ref[indel[0] + range[0] - 1] == ref[indel[1] + range[0] - 1]:
        range[0] = range[0] - 1
    while indel[0] + range[1] < len(ref) and indel[1] + range[1] < len(ref) and ref[indel[0] + range[1]] == ref[indel[1] + range[1]]:
        range[1] = range[1] + 1
    return range

def indel_dis(indel1, indel2, ref):
    range1, range2 = get_range(indel1, ref), get_range(indel2, ref)
    return NDsegment_distance(indel1, indel2, range1, range2)

ref = sys.argv[1]

for line in sys.stdin:
    mindis = numpy.inf
    query, refl1, refr1, queryl1, queryr1, c2 = line.split("\t", 5)
    indel1 = numpy.array([int(refl1), int(refr1), int(queryl1), int(queryr1)])
    c2 = c2.split("\t")
    for i in range(0, len(c2), 4):
        indel2 = numpy.array([int(c2[i]), int(c2[i + 1]), int(c2[i + 2]), int(c2[i + 3])])
        disv = indel_dis(indel1, indel2, ref)
        dis = numpy.sum(numpy.abs(disv))
        if dis < mindis:
            mindis = dis
            mindisv = disv
    sys.stdout.write(f"{query}\t{mindisv[0]}\t{mindisv[1]}\t{mindisv[2]}\t{mindisv[3]}\n")