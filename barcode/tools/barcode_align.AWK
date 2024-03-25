#!/app/gawk-5.3.0/gawk -f
# Usage: barcode_align.AWK <fqR1.demultiplex -- project_path fqR1 ext1up ext2up totalCount

function execRearr(ref1len, spliter2seqGroup)
{
    rn = 0
    while (cmd ref1len | getline result)
    {
        print result > fqR1 ".alg"
        ++rn
        if (rn % 3 != 1)
            continue
        split(result, fields, "\t")
        print spliter2seqGroup, result, ext1up, ref1len + ext2up, sprintf("%.2f%", fields[2] / totalCount * 100) > fqR1 ".table"
    }
    close(cmd ref1len)
}

BEGIN{
    FS = "\t"
    OFS = "\t"
    project_path = ARGV[1]
    fqR1 = ARGV[2]
    ext1up = ARGV[3]
    ext2up = ARGV[4]
    totalCount = ARGV[5] 
    for (i = 1; i <= 5; ++i)
        delete ARGV[i]
    
    cmd = project_path "/Rearrangement/build/rearrangement <" fqR1 ".countfile 3<" fqR1 ".reference -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | " project_path "/gawk-5.3.0/gawk -f" project_path "/tools/correct_micro_homology.AWK -- " ext1up " NGG " ext2up " NGG "
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
        {
            close(fqR1 ".countfile")
            execRearr(length(ref1), spliter2seqGroup)
            printf("\n") > fqR1 ".alg"
        }
        spliter2seqGroup = spliter2seq
        printf("0\n%s\n%d\n%d\n%s\n%d\n", ref1, ext1up, ext2up, ref2, length(ref2)) > fqR1 ".reference"
        close(fqR1 ".reference")
        print spliter2seqGroup > fqR1 ".alg"
    }
    printf("%s\t%d\n", substr(R2CutAdapt, endOfSpliter2InR2 + 4), count) > fqR1 ".countfile"
}
END{
    close(fqR1 ".countfile")
    execRearr(length(ref1), spliter2seqGroup)
    printf("\n") > fqR1 ".alg"
    system("rm " fqR1 ".countfile")
    system("rm " fqR1 ".reference")
}