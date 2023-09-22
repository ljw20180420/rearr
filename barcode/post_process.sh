#!/bin/bash
get_indel()
{
    awk -F "\t" -v OFS="\t" '
        NR==1{print $0, "left_del", "right_del", "temp_left_ins", "temp_right_ins", "random_ins", "indel_type"}
        NR>1{
            ldel = ($10 > $6 ? $10 - $6 : 0);
            rdel = ($8 > $11 ? $8 - $11 : 0);
            del = ldel + rdel;
            tlins = ($6 > $10 ? $6 - $10 : 0);
            trins = ($11 > $8 ? $11 - $8 : 0);
            rins = length($7);
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
    awk -F "\t" -v OFS="\t" -v tg=$1 'NR>1{count[$1][$tg] += $3; tgvs[$tg] = 0;}
    END{
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

# get_indel < barcode/A2_TEST.fq.final_hgsgrna_libb_all_0811-NGG.csv.barcode.table
# get_indel < barcode/A2_TEST.fq.final_hgsgrna_libb_all_0811-NGG.csv.barcode.table | get_table 17

for file in $(ls path/*.suffix)
do
    get_table col_num < $file
done
