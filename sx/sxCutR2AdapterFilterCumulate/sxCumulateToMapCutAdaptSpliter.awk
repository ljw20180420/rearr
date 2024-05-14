#!/usr/bin/env -S gawk -f

# Usage: cumulateToMapCutAdaptSplite.awk <toMapCutAdaptSpliter|count|spliter
# Cumulate the adjacent duplicate toMapCutAdaptSpliter count
# Input: toMapCutAdaptSpliter|count|spliter
# Output: toMapCutAdaptSpliter|count|spliter

BEGIN{
    FS = "\t"
    OFS = "\t"
}
{
    if ($1 != toMapCutAdaptSpliter)
    {
        if (NR > 1)
            print toMapCutAdaptSpliter, count, spliter
        toMapCutAdaptSpliter = $1
        count = $2
        spliter = $3
    }
    else
        count += $2
}
END{
    print toMapCutAdaptSpliter, count, spliter
}