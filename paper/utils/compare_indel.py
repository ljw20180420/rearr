#!/usr/bin/env python

import sys
import numpy as np


def NDsegment_correct(indel_p, indel_t, range_p, range_t):
    range = [range_p[0] - range_t[1], range_p[1] - range_t[0]]
    target = np.median(indel_t - indel_p)
    actual = (
        range[1] if target > range[1] else range[0] if target < range[0] else target
    )
    return indel_p + actual


def get_range(indel, ref1, ref2):
    range = [0, 0]
    if indel[0] > len(ref1) or indel[1] < 0:
        return range
    while (
        indel[0] + range[0] - 1 >= 0
        and indel[1] + range[0] - 1 >= 0
        and ref1[indel[0] + range[0] - 1] == ref2[indel[1] + range[0] - 1]
    ):
        range[0] = range[0] - 1
    while (
        indel[0] + range[1] < len(ref1)
        and indel[1] + range[1] < len(ref2)
        and ref1[indel[0] + range[1]] == ref2[indel[1] + range[1]]
    ):
        range[1] = range[1] + 1
    return range


def indel_correct(indel_p, indel_t, ref1, ref2):
    range_p, range_t = get_range(indel_p, ref1, ref2), get_range(indel_t, ref1, ref2)
    return NDsegment_correct(indel_p, indel_t, range_p, range_t)


ref1, ref2 = sys.argv[1], sys.argv[2]

for line in sys.stdin:
    fields = line.rstrip().split("\t")
    seq_name = fields[0]
    poss = [int(pos) for pos in fields[1:]]
    indel_t = np.array(poss[:4])
    indel_ps = np.array(poss[4:]).reshape(-1, 4)

    indel_ps_correct = np.stack(
        [indel_correct(indel_p, indel_t, ref1, ref2) for indel_p in indel_ps]
    )
    min_idx = np.abs(indel_ps_correct - indel_t).sum(axis=1).argmin()
    min_indel_p = indel_ps_correct[min_idx]

    sys.stdout.write(
        "{0}\t{1}\n".format(
            seq_name,
            "\t".join(min_indel_p.astype(str)),
        )
    )
