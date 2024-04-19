#!/usr/bin/env -S gawk -f

# Usage: endOfSpliterPos.awk <samFile
# Get the end position of spliter in reads
# Input: samFile without header
# Output: flag|spliter|endOfSpliterPos

BEGIN{
    FS = "\t"
    OFS = "\t"
}

{
    if ($6 == "*")
        print $2, $3, 0
    else
    {
        n = patsplit($6, cigarSegs, /[0-9]+[MIDNSHPX=]/)
        endOfSpliterPos = 0
        for (i = 1; i <= n; ++i)
        {
            patsplit(cigarSegs[i], num, /[0-9]+/, labels)
            if (labels[1] ~ /[MI=X]/ || labels[1] ~ /[SH]/ && i == 1)
                endOfSpliterPos += num[1]
        }
        print $2, $3, endOfSpliterPos
    }
}