#!/usr/bin/env -S gawk -f

# Usage: correct_micro_homology.AWK -- refFile NGGCCNtype1 NGGCCNtype2
# correct the microhomology to fit the 5' overhang

BEGIN{
    FS = "\t"
    refFile = ARGV[1]
    NGGCCNtype1 = ARGV[2]
    NGGCNNtype2 = ARGV[3]
    ref_id = 0
    while (getline ref < refFile)
    {
        split(ref, arr, "\t")
        cut1s[ref_id] = arr[3]
        cut2s[ref_id] = length(arr[2]) + arr[4]
        ++ref_id
    }
    for (i = 1; i <= 3; ++i)
        delete ARGV[i]
}
{
    idx = $1
    count = $2
    score = $3
    ref_id = $4
    getline refline
    getline queryline

    refstart = match(refline, /[acgtn]/) - 1
    jpos = refstart + 1 + match(substr(refline, refstart + 2), /[acgtn]/)
    jpos2 = jpos + match(substr(refline, jpos + 1), /[acgtn]/) - 1
    refend = jpos2 + 1 + match(substr(refline, jpos2 + 2), /[acgtn]/)

    uw1 = refstart
    target = substr(queryline, refstart + 1, jpos - refstart)
    gsub(/[- ]/, "", target)
    uw2 = uw1 + length(target)
    dw1 = uw2 + jpos2 - jpos
    target = substr(queryline, jpos2 + 1, refend - jpos2)
    gsub(/[- ]/, "", target)
    dw2 = dw1 + length(target)
    if (uw1 == uw2 && dw1 == dw2)
        next
    if (dw1 < dw2)
    {
        segd1 = jpos2 + match(substr(queryline, jpos2 + 1), /[ACGTN]/) - 1
        target = substr(refline, 1, segd1)
        gsub(/[- ]/, "", target)
        ds1 = length(target)
        segd2 = jpos2 + match(substr(queryline, jpos2 + 1, refend - jpos2), /[ACGTN][- ]*$/)
        target = substr(refline, segd1 + 1, segd2 - segd1)
        gsub(/[- ]/, "", target)
        ds2 = ds1 + length(target)
        if (uw1 == uw2)
        {
            us1 = us2 = (ds1 < cut2 ? (cut1 - (cut2 - ds1)) : cut1)
            segu1 = segu2 = us1 + refstart
        }
    }
    if (uw1 < uw2)
    {
        segu1 = refstart + match(substr(queryline, refstart + 1), /[ACGTN]/) - 1
        target = substr(refline, 1, segu1)
        gsub(/[- ]/, "", target)
        us1 = length(target)
        segu2 = match(substr(queryline, 1, jpos), /[ACGTN][- ]*$/)
        target = substr(refline, segu1 + 1, segu2 - segu1)
        gsub(/[- ]/, "", target)
        us2 = us1 + length(target)
        if (dw1 == dw2)
        {
            ds1 = ds2 = (us2 > cut1 ? (cut2 + (us2 - cut1)) : cut2)
            target = refline
            gsub(/[- ]/, "", target)
            reflen = length(target)
            segd1 = segd2 = refend - (reflen - ds1)
        }
    }
    if (jpos == jpos2)
    {
        split(refline, refarr, "")
        split(queryline, queryarr, "")
        if (NGGCCNtype1 == "NGG")
            itersize = (us2 > cut1 ? (us2 - cut1) : (cut1 - us2))
        else if (NGGCCNtype2 == "CCN")
            itersize = (ds1 > cut2 ? (ds1 - cut2) : (cut2 - ds1))
        else
            itersize = 0
        for (i = 0; i < itersize; ++i)
        {
            if (NGGCCNtype1 == "NGG" && us2 < cut1 || NGGCCNtype1 == "CCN" && NGGCCNtype2 == "CCN" && ds1 < cut2)
            {
                if (dw1 < dw2 && segu2 < jpos && toupper(refarr[segu2 + 1]) == toupper(refarr[segd1 + 1]))
                {
                    queryarr[segu2 + 1] = queryarr[segd1 + 1]
                    queryarr[segd1 + 1] = "-"
                    ++segu2
                    ++segd1
                    ++us2
                    ++ds1
                    ++uw2
                    ++dw1
                }
                else
                {
                    if (dw1 == dw2)
                        ds1 = ds2 = (us2 > cut1 ? (cut2 + (us2 - cut1)) : cut2)
                    break
                }
            }
            if (NGGCCNtype1 == "NGG" && us2 > cut1 || NGGCCNtype1 == "CCN" && NGGCCNtype2 == "CCN" && ds1 > cut2)
            {
                if (uw1 < uw2 && segd1 > jpos && toupper(refarr[segu2]) == toupper(refarr[segd1]))
                {
                    queryarr[segd1] = queryarr[segu2]
                    queryarr[segu2] = "-"
                    --segu2
                    --segd1
                    --us2
                    --ds1
                    --uw2
                    --dw1
                }
                else
                {
                    if (uw1 == uw2)
                        us1 = us2 = (ds1 < cut2 ? (cut1 - (cut2 - ds1)) : cut1)
                    break
                }
            }
        }
        queryline = ""
        for (s = 1; s <= length(queryarr); ++s)
            queryline = queryline queryarr[s]
    }
    printf("%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n%s\n%s\n", idx, count, score, ref_id, substr(queryline, 1, refstart), us1, uw1, us2, uw2, substr(queryline, jpos + 1, jpos2 - jpos), ds1, dw1, ds2, dw2, substr(queryline, refend + 1), cut1s[ref_id], cut2s[ref_id], refline, queryline)
}
