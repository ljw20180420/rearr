#!/usr/bin/env -S gawk -f
# Usage: barcode_align.AWK -- project_path fqR1 ext1up ext2up totalCount
# stdin: fqR1.demultiplex
# Description: align read batches over spliter2seq (reads in the sam batch have the same spliter2seq)

function execRearr(ref1len, spliter2seqGroup)
{
    # append alignment results to fqR1.alg and fqR1.table
    close(fqR1 ".countfile") # close to flush the buffer to fqR1.countfile
    cut2 = ref1len + ext2up
    rn = 0
    while (cmd ref1len | getline result)
    {
        print result > fqR1 ".alg" # append all lines to fqR1.alg
        ++rn
        if (rn % 3 != 1) # only header lines (1, 4, 7 ...) are write to fqR1.table
            continue
        split(result, fields, "\t")
        count = fields[2]
        percent = sprintf("%.2f%", fields[2] / totalCount * 100)
        print spliter2seqGroup, result, cut1, cut2, percent > fqR1 ".table"
    }
    close(cmd ref1len) # awk need close pipe before next open
    printf("\n") > fqR1 ".alg" # end the alignment batch by an empty line
}

BEGIN{
    FS = "\t"
    OFS = "\t"
    project_path = ARGV[1]
    fqR1 = ARGV[2]
    ext1up = ARGV[3]
    ext2up = ARGV[4]
    totalCount = ARGV[5]
    cut1 = ext1up
    b3t5 = 3 # distance between barcode 3' and target 5' in R2
    for (i = 1; i <= 5; ++i)
        delete ARGV[i]
    
    cmd = project_path "/Rearrangement/build/rearrangement <" fqR1 ".countfile 3<" fqR1 ".reference -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | " project_path "/tools/correct_micro_homology.AWK -- " ext1up " NGG " ext2up " NGG "
    printf("spliter2\tindex\tcount\tscore\tupdangle\tref_start1\tquery_start1\tref_end1\tquery_end1\trandom_insertion\tref_start2\tquery_start2\tref_end2\tquery_end2\tdowndangle\tcut1\tcut2\tpercent\n") > fqR1 ".table"
}

{
    count = $2
    R2CutAdapt = $3
    endOfSpliter2InR2 = $4
    spliter2seq = $5
    ref1 = $6
    ref2 = $7
    if (spliter2seqGroup != spliter2seq)
    {
        if (spliter2seqGroup)
            execRearr(length(ref1), spliter2seqGroup)
        spliter2seqGroup = spliter2seq
        printf("0\n%s\n%d\n%d\n%s\n%d\n", ref1, ext1up, ext2up, ref2, length(ref2)) > fqR1 ".reference" # it is assumed that reads have the same spliter2seq also have the same ref1, ref2, cut1, cut2
        close(fqR1 ".reference")
        print spliter2seqGroup > fqR1 ".alg" # start the alignment batch by spliter2seq
    }
    printf("%s\t%d\n", substr(R2CutAdapt, endOfSpliter2InR2 + b3t5 + 1), count) > fqR1 ".countfile"
}
END{
    execRearr(length(ref1), spliter2seqGroup)
    system("rm " fqR1 ".countfile")
    system("rm " fqR1 ".reference")
}