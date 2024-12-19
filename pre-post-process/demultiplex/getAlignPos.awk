#!/usr/bin/env -S gawk -f

# Usage: getAlignPos.awk <samFile
# Get the start\end position of local alignment for spliter\sequence
# Input: samFile without header
# Output: flag|spliter|spliterStart|spliterEnd|seqStart|seqEnd

BEGIN{
    FS = "\t"
    OFS = "\t"
}

{
    if ($6 == "*") {
        print $2, $3, 0, 0, 0, 0
    } else {
        n = patsplit($6, cigarSegs, /[0-9]+[MIDNSHPX=]/)
        spliterStart = $4 - 1 # sam file is 1-based, so minus 1
        spliterEnd = spliterStart 
        seqStart = 0
        seqEnd = 0
        for (i = 1; i <= n; ++i) {
            patsplit(cigarSegs[i], num, /[0-9]+/, labels)
            if (labels[1] ~ /[MI=X]/ || labels[1] ~ /[SH]/ && i == 1) {
                seqEnd += num[1]
                if (labels[1] ~ /[SH]/ && i == 1) {
                    seqStart += num[1]
                }
            }
            if (labels[1] ~ /[MD=X]/) {
                spliterEnd += num[1]
            }
        }
        print $2, $3, spliterStart, spliterEnd, seqStart, seqEnd
    }
}