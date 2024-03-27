#!/usr/bin/env -S gawk -f
# Usage: align_align.AWK -- ref1len cut1 cut2
# stdin: fqR1.alg
# align the cut points and random insertion of alignments
# cpos1/2 and jpos are gapped positions (including random insertion gap of length midlen)
# there are exactly cut1/2 bases upstream to cpos1/2
# cpos1 is pushed as upstream as possible; cpos2 is pushed as downstream as possible
# nNNNN-NNNNNNNN|||||NNNN-NNNNn||||-------nNNNN-NNNN|||||NNNN-NNNNNNNNn
#               cpos1          jpos midlen          cpos2


function repeatChar(repchar, repn)
{
    if (repn == 0)
        return ""
    if (repn == 1)
        return repchar
    log2repn = int(log(repn - 1)/log(2)) + 1
    repstr = repchar
    for (lr = 1; lr <= log2repn; lr++) {
        repstr = repstr repstr
    }
    return substr(repstr, 1, repn)
}

function printANSI(queryline, jpos, midlen, cpos1, cpos2)
{
    ext1 = jpos - cpos1
    ext2 = cpos2 - jpos - midlen
    print repeatChar(" ", cpos1max - cpos1) substr(queryline, 1, cpos1) "\033[4;38;2;255;0;0m" substr(queryline, cpos1 + 1, jpos - cpos1) repeatChar(" ", ext1max - ext1) "\033[38;2;0;255;0m" substr(queryline, jpos + 1, midlen) repeatChar(" ", midmax - midlen + ext2max - ext2) "\033[4;38;2;0;0;255m" substr(queryline, jpos + midlen + 1, cpos2 - jpos -midlen) "\033[0m"  substr(queryline, cpos2 + 1)
}

BEGIN{
    FS = "\t"
    ref1len = ARGV[1]
    cut1 = ARGV[2]
    cut2 = ARGV[3]
    for (i = 1; i <= 3; ++i)
        delete ARGV[i]
}

{
    idx = int((NR - 1) / 3) + 1
}
# header line
NR % 3 == 1{
    midlens[idx] = length($9)
}
# the first reference line
NR == 2{
    reflinefirst = $0
}
# reference line
NR % 3 == 2{
    jpos = match($0, /[acgtn]/)
    jpos += match(substr($0, jpos + 1), /[acgtn]/)
    jposs[idx] = jpos

    cpos1 = cut1
    sstr = substr($0, 1, cpos1)
    gsub(/[- ]/, "", sstr)
    for (cnum = length(sstr); cnum < cut1; ++cnum)
        cpos1 += match(substr($0, cpos1 + 1), /[acgtnACGTN]/)
    cpos1s[idx] = cpos1

    cpos2 = jposs[idx] + midlens[idx] + cut2 - ref1len
    sstr = substr($0, 1, cpos2)
    gsub(/[- ]/, "", sstr)
    for (cnum = length(sstr); cnum < cut2 + 1; ++cnum)
        cpos2 += match(substr($0, cpos2 + 1), /[acgtnACGTN]/)
    cpos2s[idx] = --cpos2 # cnum = cut2 + 1, so cpos2 need to decrease 1
}
# query line
NR % 3 == 0{
    querylines[idx] = $0
}
END{
    for (i = 1; i <= idx; ++i)
    {
        cpos1max = cpos1s[i] > cpos1max ? cpos1s[i] : cpos1max
        ext1max = (jposs[i] - cpos1s[i]) > ext1max ? (jposs[i] - cpos1s[i]) : ext1max
        midmax = midlens[i] > midmax ? midlens[i] : midmax
        ext2max = (cpos2s[i] - jposs[i] - midlens[i]) > ext2max ? (cpos2s[i] - jposs[i] - midlens[i]) : ext2max
    }

    printANSI(reflinefirst, jposs[1], midlens[1], cpos1s[1], cpos2s[1])
    for (i = 1; i <= idx; ++i)
        printANSI(querylines[i], jposs[i], midlens[i], cpos1s[i], cpos2s[i])
}