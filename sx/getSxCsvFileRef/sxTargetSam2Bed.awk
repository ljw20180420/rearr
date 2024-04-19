#!/usr/bin/env -S gawk -f

# Usage: sxTargetSam2Bed.awk <targetSamFile -- ext1up ext1down ext2up ext2down
# get the cut position according to the target mapping (by bowtie2), and expand the cut position to genome range of ref1 and ref2 according to ext parameters

BEGIN{
    FS = "\t"
    OFS = "\t"
    ext1up = ARGV[1]
    ext1down = ARGV[2]
    ext2up = ARGV[3]
    ext2down = ARGV[4]
    for (i = 1; i <= 4; ++i)
        delete ARGV[i]
}

{
    qname = $1
    flag = $2
    if (int(flag / 16) % 2)
        refstrand = "+"
    else
        refstrand = "-"
    chr = $3
    pos = $4 - 1
    if (refstrand == "+")
    {
        target_length = substr($6, 1, length($6) - 1)
        cut = pos + target_length - 16 - 6
        print chr, cut - ext1up, cut + ext1down, qname, ".", "+"
        print chr, cut - ext2up, cut + ext2down, qname, ".", "+"
    }
    else
    {
        cut = pos + 16 + 6
        print chr, cut - ext1down, cut + ext1up, qname, ".", "-"
        print chr, cut - ext2down, cut + ext2up, qname, ".", "-"
    }
}