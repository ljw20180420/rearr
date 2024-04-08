#!/usr/bin/env -S gawk -f
# Usage: align_align.AWK -- ref1len cut1 cut2
# stdin: fqR1.alg
# align the cut points and random insertion of alignments
# cpos1/2 and jpos are gapped positions (including random insertion gap of length midlen)
# there are exactly cut1/2 bases upstream to cpos1/2
# cpos1 is pushed as upstream as possible; cpos2 is pushed as downstream as possible
# nNNNNNNNNNNNN|NNNNNNNNn-------nNNNNNNNN|NNNNNNNNNNNNn
#            cpos1     jpos  +midlen   cpos2

function getGapAnsiCodes(i)
{
    CrossCut1[idx] = 0
    CrossCut2[idx] = 0
    for (i = 1; i <= mmns[idx]; ++i)
    {
        if (i == 1 || gapposss[idx][i - 1] == ref1len || gapposss[idx][i] - gapposss[idx][i - 1] > 1 || gapposss[idx][i] == ref1len || i == mmns[idx])
        {
            if (i > 1 && gapposss[idx][i - 1] != ref1len)
                GapAnsiCodess[idx][gapposss[idx][i - 1]] = "\033[49m"
            if (gapposss[idx][i] != ref1len && i != mmns[idx])
                GapAnsiCodess[idx][gapposss[idx][i] - 1] = "\033[46m"
        }
    }
}

function modifyGapAnsiCodes(GapAnsiCodes, key)
{
    delete GapAnsiCodesPrint
    for (key in GapAnsiCodes)
    {
        keyPrint = (key >= ref1len) ? (key + midMax) : key 
        GapAnsiCodesPrint[keyPrint] = GapAnsiCodes[key]
    }
}

function getCutJuncAnsiCodes()
{
    CutJuncAnsiCodes[cut1] = "\033[31m"
    CutJuncAnsiCodes[ref1len] = "\033[39m\033[32m"
    CutJuncAnsiCodes[ref1len + midMax] = "\033[39m\033[33m"
    CutJuncAnsiCodes[cut2 + midMax] = "\033[39m"
}

function getInsertions(refline, queryline, i)
{
    mmns[idx] = split(refline, mms, /-+/ , gapss[idx])
    gapposss[idx][0] = 0
    for (i = 1; i <= mmns[idx]; ++i)
        gapposss[idx][i] = gapposss[idx][i - 1] + length(mms[i])
    querylineRmIns = ""
    mmend = 0
    for (i = 1; i <= mmns[idx]; ++i)
    {
        mmstart = mmend + length(gapss[idx][i - 1])
        mmend = mmstart + length(mms[i])
        querylineRmIns = querylineRmIns substr(queryline, mmstart + 1, mmend - mmstart)
    }
    return querylineRmIns
}

function repeatChar(repchar, repn, i)
{
    if (repn == 0)
        return ""
    if (repn == 1)
        return repchar
    log2repn = int(log(repn - 1)/log(2)) + 1
    repstr = repchar
    for (i = 1; i <= log2repn; i++) {
        repstr = repstr repstr
    }
    return substr(repstr, 1, repn)
}

function printANSI(queryline, gapposs, mmn, mid, N, n, i)
{
    queryline = substr(queryline, 1, ref1len) mid repeatChar(" ", midMax - length(mid), 0) substr(queryline, ref1len + 1)
    N = split(queryline, querylineArr, "")
    outputline = ""
    i = 1
    for (n = 1; n <= N; ++n)
    {
        if (i < mmn && n == gapposs[i] && n != ref1len)
        {
            outputline = outputline "(" querylineArr[n] ")"
            ++i
        }
        else
        {
            outputline = outputline "{" querylineArr[n] "}"
        }
        if (n == cut1)
            outputline = outputline "\033[31m"
        else if (n == ref1len)
            outputline = outputline "\033[39m\033[32m"
        else if (n == ref1len + midMax)
            outputline = outputline "\033[39m\033[33m"
        else if (n == cut2 + midMax)
            outputline = outputline "\033[39m"
    }
    print outputline
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
    mids[idx] = $9
}
# the first reference line
NR == 2{
    ref = $0
    gsub(/[- ]/, "", ref)
}
# reference line
NR % 3 == 2{
    refline = $0
}
# query line
NR % 3 == 0{
    querylines[idx] = getInsertions(refline, $0, 0)
    getGapAnsiCodes(0)
}
END{
    for (i = 1; i <= idx; ++i)
        midMax = length(mids[i]) > midMax ? length(mids[i]) : midMax
    getCutJuncAnsiCodes()
    
    split("", emptyarray)
    split("", GapAnsiCodesPrint)
    printANSI(ref, emptyarray, 0, repeatChar(" ", midMax, 0), 0, 0, 0)
    for (i = 1; i <= idx; ++i)
    {
        modifyGapAnsiCodes(GapAnsiCodess[i], 0)
        printANSI(querylines[i], gapposss[i], mmns[i], mids[i], 0, 0, 0)
    }
}