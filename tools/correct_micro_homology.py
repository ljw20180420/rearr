#!/usr/bin/env python
import sys, re, more_itertools
def correct_micro(cut1, cut2, NGGCCNtype1, NGGCCNtype2, fd):
    # modify the alignment within micro homology such that the upstream end is blunt for NGG and downstream end is blunt for CCN if possible
    relow = re.compile("[acgtn]")
    reup = re.compile("[ACGTN]")
    for header, refline, queryline in more_itertools.batched(fd, 3):
        idx, count, score = header.rstrip("\n").split("\t", 3)
        refline, queryline = refline.rstrip(), queryline.rstrip()
        refstart = relow.search(refline, pos = 0).span()[0]
        jpos = relow.search(refline, pos = refstart + 1).span()[1]
        jpos2 = relow.search(refline, pos = jpos).span()[0]
        refend = relow.search(refline, pos = jpos2 + 1).span()[1]
        uw1 = refstart
        uw2 = uw1 + len(queryline[refstart:jpos].replace(" ", "").replace("-", ""))
        dw1 = uw2 + jpos2 - jpos
        dw2 = dw1 + len(queryline[jpos2:refend].replace(" ", "").replace("-", ""))
        if uw1 == uw2 and dw1 == dw2:
            continue
        if dw1 < dw2:
            segd1 = reup.search(queryline, pos = jpos2).span()[0]
            ds1 = len(refline[:segd1].replace(" ", "").replace("-", ""))
            segd2 = refend - reup.search(queryline[refend-1:jpos2-1:-1]).span()[0]
            ds2 = ds1 + len(refline[segd1:segd2].replace(" ", "").replace("-", ""))
            if uw1 == uw2:
                us1 = us2 = cut1 - (cut2 - ds1) if ds1 < cut2 else cut1
                segu1 = segu2 = us1 + refstart
        if uw1 < uw2:
            segu1 = reup.search(queryline, pos = refstart).span()[0]
            us1 = len(refline[:segu1].replace(" ", "").replace("-", ""))
            segu2 = jpos - reup.search(queryline[jpos-1::-1]).span()[0]
            us2 = us1 + len(refline[segu1:segu2].replace(" ", "").replace("-", ""))
            if dw1 == dw2:
                ds1 = ds2 = cut2 + (us2 - cut1) if us2 > cut1 else cut2
                reflen = len(refline.replace(" ", "").replace("-", ""))
                segd1 = segd2 = refend - (reflen - ds1)
        if jpos == jpos2:
            queryline = list(queryline)
            itersize = abs(us2 - cut1) if NGGCCNtype1 == "NGG" else abs(ds1 -cut2) if NGGCCNtype2 == "CCN" else 0
            for _ in range(itersize):
                if NGGCCNtype1 == "NGG" and us2 < cut1 or NGGCCNtype1 == "CCN" and NGGCCNtype2 == "CCN" and ds1 < cut2:
                    if dw1 < dw2 and segu2 < jpos and refline[segu2].upper() == refline[segd1].upper():
                        queryline[segu2] = queryline[segd1]
                        queryline[segd1] = "-"
                        segu2, segd1, us2, ds1, uw2, dw1 = segu2 + 1, segd1 + 1, us2 + 1, ds1 + 1, uw2 + 1, dw1 + 1
                    else:
                        if dw1 == dw2:
                            ds1 = ds2 = cut2 + (us2 - cut1) if us2 > cut1 else cut2
                        break
                if NGGCCNtype1 == "NGG" and us2 > cut1 or NGGCCNtype1 == "CCN" and NGGCCNtype2 == "CCN" and ds1 > cut2:
                    if uw1 < uw2 and segd1 > jpos and refline[segu2 - 1].upper() == refline[segd1 - 1].upper():
                        queryline[segd1 - 1] = queryline[segu2 - 1]
                        queryline[segu2 - 1] = "-"
                        segu2, segd1, us2, ds1, uw2, dw1 = segu2 - 1, segd1 - 1, us2 - 1, ds1 - 1, uw2 - 1, dw1 - 1
                    else:
                        if uw1 == uw2:
                            us1 = us2 = cut1 - (cut2 - ds1) if ds1 < cut2 else cut1
                        break
            queryline = "".join(queryline)
        yield f"{idx}\t{count}\t{score}\t{queryline[:refstart]}\t{us1}\t{uw1}\t{us2}\t{uw2}\t{queryline[jpos:jpos2]}\t{ds1}\t{dw1}\t{ds2}\t{dw2}\t{queryline[refend:]}", refline, queryline
        

if __name__ == "__main__":
    cut1, NGGCCNtype1, cut2, NGGCCNtype2, ref1len, = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
    for header, refline, queryline in correct_micro(cut1, ref1len + cut2, NGGCCNtype1, NGGCCNtype2, sys.stdin):
        sys.stdout.write(f"{header}\n")
        sys.stdout.write(f"{refline}\n")
        sys.stdout.write(f"{queryline}\n")
