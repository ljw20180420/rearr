#!/bin/bash
# Usage
# Get a long table of indel
# rearr_barcode_post_process.sh <fqR1.table
# To hist a certain column
# rearr_barcode_post_process.sh <fqR1.table column

get_indel()
{
    gawk -F "\t" -v OFS="\t" '
        NR==1{
                print $0, "left_del", "right_del", "temp_left_ins", "temp_right_ins", "random_ins", "indel_type"
        }
        NR>1{
            ref_end1 = $8
            random_insertion = $10
            ref_start2 = $11
            cut1 = $16
            cut2 = $17

            left_del = (cut1 > ref_end1 ? cut1 - ref_end1 : 0)
            right_del = (ref_start2 > cut2 ? ref_start2 - cut2 : 0)
            del = left_del + right_del

            temp_left_ins = (ref_end1 > cut1 ? ref_end1 - cut1 : 0)
            temp_right_ins = (cut2 > ref_start2 ? cut2 - ref_start2 : 0)
            random_ins = length(random_insertion)
            ins = temp_left_ins + temp_right_ins + random_ins

            if (del > 0 && ins > 0) 
                indel_type = "indel"
            else if (del > 0 && ins == 0)
                indel_type = "del"
            else if (del == 0 && ins > 0)
                indel_type = "ins"
            else
                indel_type = "WT"

            print $0, left_del, right_del, temp_left_ins, temp_right_ins, random_ins, indel_type
        }
    '
}

get_table()
{
    gawk -v FS="\t" -v tg=$1 'NR>1{count[$1][$tg] += $3; tgvs[$tg] = 0;}
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

project_path="$(dirname $(realpath $0))/../.."

if [ -n "$1" ]
then
    get_indel | get_table $1
else
    get_indel
fi
