#!/bin/bash
get_indel()
{
    awk -F "\t" -v OFS="\t" '
        {
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

input=${1:-"ljhlyz/AN1-SG4-M1B-1-1_R1.fq.gz"} # the input file
# case $input in # count duplicate reads, support fasta, fastq, and their compressions
#     *.fq.gz|*.fastq.gz)
#         echo "the input is gzip fastq file"
#         zcat $input | sed -n '2~4p' | sort | uniq -c | sort -k1,1nr | awk -v OFS="\t" '{print $2, $1}' >$input.count;;
#     *.fa.gz|*.fasta.gz)
#         echo "the input is gzip fasta file"
#         zcat $input | sed -n '2~2p' | sort | uniq -c | sort -k1,1nr | awk -v OFS="\t" '{print $2, $1}' >$input.count;;
#     *.fq|*.fastq)
#         echo "the input is fastq file"
#         sed -n '2~4p' | sort | uniq -c | sort -k1,1nr | awk -v OFS="\t" '{print $2, $1}' >$input.count;;
#     *.fa|*.fasta)
#         echo "the input is fasta file"
#         sed -n '2~2p' | sort | uniq -c | sort -k1,1nr | awk -v OFS="\t" '{print $2, $1}' >$input.count;;
# esac
ref=${2:-"CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGTTGCTGTTGCTGGTGCTGATGGTGATGTGTTGAGACTGGTGGGTGGGCGGTGGACTGGGCCCCAGTAGAGGGAGGGAAGGGGCCTGGATGGGCATTGCTGTT"} # reference
sgRNA=${3:-"GGTGATGTGTTGAGACTGGT"} # sgRNA
exec=${4:-"Rearrangement/build/rearrangement"} # executable for alignment
ext1=${5:-30} # upstream end downstream extension for template inserion (default: 30)
ext2=${6:-30} # downstream end upstream extension (default: 30)

cut=$(generate_ref_file.py $input.count $ref $sgRNA $exec $ext1 $ext2) # prepare the reference file and return the cut point

$exec -file $input.count -ref_file ref_file -ALIGN_MAX 1 -THR_NUM 24 -u -1 -v -3 -s0 -2 -qv -3 | sed -nr 'N;N;s/\n/\t/g;p' | sort -k1,1n | awk -F "\t" '{for (i=1; i<=NF-3; ++i) printf("%s\t",$i); printf("%s\n%s\n%s\n", $(NF-2), $(NF-1), $NF);}' | tee $input.alg.$cut.$ext1.$ext2 | correct_micro_homology.py $(($cut + $ext1)) $cut $(($cut + $ext1 + $ext2)) | tee $input.correct.$cut.$ext1.$ext2 | awk -v OFS="\t" -v cut1=$cut -v cut2=$(($cut + $ext1 + $ext2)) 'NR%3==1{print $0, cut1, cut2}' | get_indel > $input.table.$cut.$ext1.$ext2 # align reads (input.alg), correct micro homology (input.correct)