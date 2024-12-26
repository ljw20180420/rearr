#!/usr/bin/env -S gawk -f

# Usage: sxCumulateToMapCutAdaptSplite.awk <query|count|refId
# Accumulate the counts of adjacent duplicated queries
# Input: query|count|refId
# Output: query|count|refId

BEGIN{
    FS = "\t"
    OFS = "\t"
}
{
    if ($1 != toMapCutAdaptSpliter)
    {
        if (NR > 1)
            print toMapCutAdaptSpliter, count, refId
        toMapCutAdaptSpliter = $1
        count = $2
        refId = $3
    }
    else
        count += $2
}
END{
    print toMapCutAdaptSpliter, count, refId
}
