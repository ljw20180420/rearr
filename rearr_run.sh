#!/bin/bash
# usage: rearr_run.sh input ref1 ref2 cut1 cut2 NGGCCNtype1 NGGCCNtype2
# input can be (compressed) fasta/fastq

project_path="$(dirname $(realpath $0))"

input=$1 # the input file
case $input in # count duplicate reads, support fasta, fastq, and their compressions
    *.fq.gz|*.fastq.gz)
        echo "the input is gzip fastq file"
        $project_path/pv-1.8.5/pv -c -N "count $input" "$input" | zcat | sed -n '2~4p' | sort | uniq -c | sort -k1,1nr | gawk -v OFS="\t" '{print $2, $1}' >"$input.count";;
    *.fa.gz|*.fasta.gz)
        echo "the input is gzip fasta file"
        $project_path/pv-1.8.5/pv -c -N "count $input" "$input" | zcat | sed -n '2~2p' | sort | uniq -c | sort -k1,1nr | gawk -v OFS="\t" '{print $2, $1}' >"$input.count";;
    *.fq|*.fastq)
        echo "the input is fastq file"
        $project_path/pv-1.8.5/pv -c -N "count $input" "$input" | sed -n '2~4p' | sort | uniq -c | sort -k1,1nr | gawk -v OFS="\t" '{print $2, $1}' >"$input.count";;
    *.fa|*.fasta)
        echo "the input is fasta file"
        $project_path/pv-1.8.5/pv -c -N "count $input" "$input" | sed -n '2~2p' | sort | uniq -c | sort -k1,1nr | gawk -v OFS="\t" '{print $2, $1}' >"$input.count";;
esac
ref1=$2 # reference1
ref2=$3 # reference2
cut1=$4 # cut point of reference1
cut2=$5 # cut point of reference2
NGGCCNtype1=$6 # NGGCCNtype of reference1
NGGCCNtype2=$7 # NGGCCNtype of reference2

printf "0\n%s\n%d\n%d\n%s\n%d\n" "$ref1" "$cut1" "$cut2" "$ref2" ${#ref2} >"$input.ref.$cut1.$cut2.${#ref1}"

$project_path/pv-1.8.5/pv -c -N "align $input" "$input.count" | $project_path/Rearrangement/build/rearrangement 3<"$input.ref.$cut1.$cut2.${#ref1}" -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | $project_path/tools/correct_micro_homology.AWK -- "$cut1" "$NGGCCNtype1" "$cut2" "$NGGCCNtype2" ${#ref1} | tee "$input.alg.$cut1.$cut2.${#ref1}" | gawk -F "\t" -v OFS="\t" -v cut1="$cut1" -v cut2=$((${#ref1} + "$cut2")) -v total_count="$(gawk -F "\t" '{count += $2} END{print count}' $input.count)" '
    BEGIN{
        print "index", "count", "score", "updangle", "ref_start1", "query_start1", "ref_end1", "query_end1", "random_insertion", "ref_start2","query_start2", "ref_end2", "query_end2", "downdangle", "cut1", "cut2", "percent"
    }
    NR%3==1{
        printf("%s\t%d\t%d\t%.2f%\n", $0, cut1, cut2, $2/total_count*100)
        }
' >"$input.table.$cut1.$cut2.${#ref1}" # align reads (input.alg), correct micro homology (input.correct)