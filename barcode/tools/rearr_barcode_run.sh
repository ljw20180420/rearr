#!/bin/bash

# Usage: printf "%s\n" fastq1 fastq2 ... | rearr_barcode_run.sh

find_barcode()
{
    local minscoreR1=100
    local minscoreR2=30
    local fq1=$1
    local csvfile=$2
    local fq2=${fq1%.*}.R2.${fq1##*.}

    pv-1.8.5/pv -c -N "count $fq1" "$fq2" | sed -n '2~4p' | paste - <(sed -n '2~4p' "$fq1") | sort | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' >"$fq1.count"

    pv-1.8.5/pv -c -N "demultiplex $fq1" "$fq1.count" | cut -f1 | bowtie2 --quiet --norc --mm --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,$minscoreR2 -r -x "$csvfile.primer+barcode" -U - 2>/dev/null | samtools view | cut -f2,3,6 | $python_exec -c '
import sys, pysam
for line in sys.stdin:
    flag, barcode, CIGAR = line.split("\t")
    pAS = pysam.AlignedSegment()
    pAS.cigarstring = CIGAR
    sys.stdout.write(f"{flag}\t{barcode}\t{pAS.query_alignment_end}\n")
    ' | awk -F "\t" -v csvfile="$csvfile" -v OFS="\t" '
        BEGIN{
            pbfa = csvfile ".primer+barcode.fa"; 
            ref12 = csvfile ".ref12"
            while (getline bn <pbfa) 
            {
                getline bc <pbfa;
                bc = substr(bc, 22)
                bn2bc[substr(bn, 2)] = bc;
                getline r12 <ref12;
                bn2r12[substr(bn, 2)] = r12
            }
        }
        {
            print $1, $2, $3, bn2bc[$2], bn2r12[$2]
        }
    ' | paste "$fq1.count" <(cut -f1 "$fq1.count" | sed '=' | sed '1~2s/^/>s/' | cutadapt -a GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC - 2> /dev/null | sed '1~2d') - <(cut -f2 "$fq1.count" | bowtie2 --quiet --norc --mm --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,$minscoreR1 -r -x "$csvfile.sgRNA+scaffold" -U - 2>/dev/null | samtools view | cut -f2,3) | awk -F "\t" -v OFS="\t" -v minscoreR2=$minscoreR2 -v fq1="$fq1" '
    {
        if (($5/4)%2 == 1 || ($11/4)%2 == 1 || $6 != $12 || $7 + minscoreR2 > length($4))
            print $0 > fq1".not_find";
        else
            print $1, $3, $4, $7, $8, $9, $10;
    }
    ' | awk -F "\t" -v OFS="\t" '
    {
        if ($1 != read2)
        {
            if (NR > 1)
                print read2, count, read2rmadapter, query_start, barcode, ref1, ref2;
            read2 = $1;
            count = $2;
            read2rmadapter = $3;
            query_start = $4;
            barcode = $5;
            ref1 = $6;
            ref2 = $7
        }
        else
            count += $2;
    }
    END{
        print read2, count, read2rmadapter, query_start, barcode, ref1, ref2;
    }
    ' | sort -k5,5 -k2,2nr
}

ext1up=50
ext2up=10
project_path="$(dirname $(realpath $0))/../.."
python_exec="$(find $project_path/renv -type f,l -name python)"

while read fq1
do
    csvfile="$($project_path/barcode/tools/infer_csvfile.sh $fq1)"
    (
        find_barcode "$fq1" "$csvfile" >"$fq1.barcode"
        pv-1.8.5/pv -c -N "align $fq1" "$fq1.barcode" | $python_exec $project_path/barcode/tools/rearr_barcode_align.py "$fq1" "$ext1up" "$ext2up"
        pv-1.8.5/pv -c -N "cal percent $fq1" "$fq1.table" | awk -F "\t" -v OFS="\t" -v total="$(tail -n+2 $fq1.table | cut -f3 | awk '{total += $0} END{print total}')" 'NR == 1{print $0, "percent"} NR > 1{printf("%s\t%.2f%\n", $0, $3 / total * 100)}' >"$fq1.table2"
        mv "$fq1.table2" "$fq1.table"
    ) &
done

wait





