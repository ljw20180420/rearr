#!/usr/bin/env python
import sys, re, more_itertools
ref1len, cut1, cut2 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
relow = re.compile("[acgtn]")
reup = re.compile("[ACGTN]")
reflines, querylines, jposs, midlens, cpos1s, cpos2s = [], [], [], [], [], []
# with open("ljhlyz/AN1-SG4-M1B-1-1_R1.fq.gz.correct", "r") as fd:
#     for header, refline, queryline in more_itertools.batched(fd, 3):
for header, refline, queryline in more_itertools.batched(sys.stdin, 3):
    reflines.append(refline)
    idx, count, score, us, ue, mid, ds, de= header.split("\t")
    midlens.append(len(mid))
    us, ue, ds, de = int(us), int(ue), int(ds), int(de)
    jpos = relow.search(refline, pos = ref1len - 1).span()[1]
    jposs.append(jpos)
    left1 = ref1len - cut1
    for cpos1 in range(jpos - 1, 0, -1):
        if refline[cpos1] != "-":
            left1 -= 1
        if left1 == 0:
            break
    while refline[cpos1 - 1] == "-":
        cpos1 -= 1
    cpos1s.append(cpos1)
    left2 = cut2 - ref1len
    for cpos2 in range(jpos + len(mid), len(refline)):
        if refline[cpos2] != "-":
            left2 -= 1
        if left2 == 0:
            break
    while refline[cpos2 + 1] == "-":
        cpos2 += 1
    cpos2 += 1
    cpos2s.append(cpos2)
    querylines.append(queryline)

cpos1max = max(cpos1s)
ext1max = max([jpos - cpos1 for jpos, cpos1 in zip(jposs, cpos1s)])
midmax = max(midlens)
ext2max = max([cpos2 - jpos - midlen for jpos, cpos2, midlen in zip(jposs, cpos2s, midlens)])

for queryline, jpos, midlen, cpos1, cpos2 in zip([reflines[0]] + querylines, [jposs[0]] + jposs, [midlens[0]] + midlens, [cpos1s[0]] + cpos1s, [cpos2s[0]] + cpos2s):
    ext1, ext2 = jpos - cpos1, cpos2 - jpos - midlen
    sys.stdout.write(" " * (cpos1max - cpos1) + queryline[:cpos1] + "\033[4;38;2;255;0;0m" + queryline[cpos1:jpos] + " " * (ext1max - ext1) + "\033[38;2;0;255;0m" + queryline[jpos:jpos + midlen] + " " * (midmax - midlen + ext2max - ext2) + "\033[4;38;2;0;0;255m" + queryline[jpos + midlen:cpos2] + "\033[0m" + queryline[cpos2:])