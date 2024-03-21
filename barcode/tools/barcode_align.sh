#!/bin/bash
# Usage: barcode_align.sh fqR1 ext1up ext2up totalCount

awk -F "\t" -v OFS="\t" -v project_path="$(dirname $(realpath $0))/../.." -v fqR1="$1" -v ext1up="$2" -v ext2up="$3" -v totalCount=$4 '
function execRearr(ref1len, spliter2seqGroup)
{
    rn = 0
    while (cmd ref1len | getline result)
    {
        print result
        ++rn
        if (rn % 3 != 1)
            continue
        split(result, fields, "\t")
        print spliter2seqGroup, result, ext1up, ref1len + ext2up, sprintf("%.2f%", fields[2] / totalCount * 100) > "/dev/fd/3"
    }
    close(cmd ref1len)
}

BEGIN{
    cmd = project_path "/Rearrangement/build/rearrangement <" fqR1 ".countfile 3<" fqR1 ".reference -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | " project_path "/tools/correct_micro_homology.AWK -- " ext1up " NGG " ext2up " NGG "
    printf("spliter2\tindex\tcount\tscore\tupdangle\tref_start1\tquery_start1\tref_end1\tquery_end1\trandom_insertion\tref_start2\tquery_start2\tref_end2\tquery_end2\tdowndangle\tcut1\tcut2\tpercent\n") > "/dev/fd/3"
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
            printf("\n")
        }
        spliter2seqGroup = spliter2seq
        printf("0\n%s\n%d\n%d\n%s\n%d\n", ref1, ext1up, ext2up, ref2, length(ref2)) > fqR1 ".reference"
        close(fqR1 ".reference")
        print spliter2seqGroup
    }
    printf("%s\t%d\n", substr(R2CutAdapt, endOfSpliter2InR2 + 4), count) > fqR1 ".countfile"
}
END{
    close(fqR1 ".countfile")
    execRearr(length(ref1), spliter2seqGroup)
    printf("\n")
    system("rm " fqR1 ".countfile")
    system("rm " fqR1 ".reference")
}
'