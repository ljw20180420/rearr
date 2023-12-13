#!/bin/bash
get_indel()
{
    awk -F "\t" -v OFS="\t" '
        NR==1{print $0, "left_del", "right_del", "temp_left_ins", "temp_right_ins", "random_ins", "indel_type"}
        NR>1{
            ldel = ($16 > $8 ? $16 - $8 : 0);
            rdel = ($11 > $17 ? $11 - $17 : 0);
            del = ldel + rdel;
            tlins = ($8 > $16 ? $8 - $16 : 0);
            trins = ($17 > $11 ? $17 - $11 : 0);
            rins = length($10);
            ins = tlins + trins + rins
            printf("%s\t%d\t%d\t%d\t%d\t%d\t", $0, ldel, rdel, tlins, trins, rins)
            if (del > 0 && ins > 0) print "indel";
            else if (del > 0 && ins == 0) print "del";
            else if (del == 0 && ins > 0) print "ins";
            else print "WT";
        }'
}

get_table()
{
    awk -v FS="\t" -v tg=$1 'NR>1{count[$1][$tg] += $3; tgvs[$tg] = 0;}
    END{
        PROCINFO["sorted_in"] = "@ind_str_asc";
        printf("barcode");
        for(kt in tgvs) 
            printf("\t%s", kt);
        printf("\n");
        for (kb in count)
        {
            printf(kb);
            for (kt in tgvs)
                printf("\t%d", count[kb][kt])
            printf("\n");
        }
    }'
}

if [ -n "$1" ]
then
    get_indel | get_table $1
else
    get_indel
fi
