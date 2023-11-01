#!/usr/bin/env python
import sys, re, more_itertools
def correct_micro(ref1len, cut1, cut2, fd):
    relow = re.compile("[acgtn]")
    reup = re.compile("[ACGTN]")
    for header, refline, queryline in more_itertools.batched(fd, 3):
        idx, count, score, _, us, ue, mid, ds, de, _= header.split("\t", 9)
        us, ue, ds, de = int(us), int(ue), int(ds), int(de)
        if us == ue and ds == de:
            continue
        if not mid and ue != cut1:
            jpos = relow.search(refline, pos = ref1len - 1).span()[1]
            if ds < de:
                qds = reup.search(queryline, pos = jpos).span()[0]
                if us == ue:
                    if ds < cut2:
                        que = us = ue = cut1 - (cut2 - ds)
                    else:
                        que = us = ue = cut1
            if us < ue:
                que = jpos - reup.search(queryline[jpos-1::-1]).span()[0]
                if ds == de:
                    if ue > cut1:
                        qds = ds = de = cut2 + (ue - cut1)
                    else:
                        qds = ds = de = cut2
            queryline = list(queryline)
            for _ in range(abs(ue - cut1)):
                if ue < cut1:
                    if refline[que].upper() == refline[qds].upper():
                        queryline[que] = queryline[qds]
                        queryline[qds] = "-"
                        que += 1
                        qds += 1
                        ue += 1
                        ds += 1
                    else:
                        break
                if ue > cut1:
                    if refline[que - 1].upper() == refline[qds - 1].upper():
                        queryline[qds - 1] = queryline[que - 1]
                        queryline[que - 1] = "-"
                        que -= 1
                        qds -= 1
                        ue -= 1
                        ds -= 1
            queryline = "".join(queryline)
        yield f"{idx}\t{count}\t{score}\t{us}\t{ue}\t{mid}\t{ds}\t{de}\n", refline, queryline
        

if __name__ == "__main__":
    ref1len, cut1, cut2 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    for header, refline, queryline in correct_micro(ref1len, cut1, cut2, sys.stdin):
        sys.stdout.write(header)
        sys.stdout.write(refline)
        sys.stdout.write(queryline)

