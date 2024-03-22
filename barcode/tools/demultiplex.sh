#!/bin/bash
# Usage: demultiplex.sh fqR1 fqR2 spliter1 spliter2 sgRNAfile ref12 minscoreR1 minscoreR2 minQueryR2
# spliter1/2 is the library to split fqR1/2
# minscoreR1/2 thres the match between fqR1/2 and spliter1/2
# ref12 is a file with side-by-side ref1 and ref2
# Remove spliter2 from 5' and adapter from 3' of R2. The remain must contain at least minQueryR2 bases.

cutadaptPlain()
{
    # Usage: cutadaptPlain <plainseq 3'adapter
    # cutadapt does not accept plainseq. This function transform plainseq to fasta before feed to cutadapt, and then transform the fasta output back to plainseq
    # Input: R2reads
    # Output: 3' trimmed R2reads
    sed '=' | sed '1~2s/^/>s/' | cutadapt -a $1 - 2> /dev/null | sed '1~2d'
}

R2map()
{
    # Usage: R2map <R2reads
    # Map R2 to spliter2, and retrieve flag, matched spliter2, and the end position of spliter2 in R2
    # Input: R2reads
    # Output: flag|spliter2|endOfSpliter2InR2
    bowtie2 --quiet --norc --mm --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,$minscoreR2 -r -x "$spliter2" -U - 2>/dev/null | samtools view | awk -F "\t" -v OFS="\t" '
    {
        if ($6 == "*")
            print $2, $3, 0
        else
        {
            n = patsplit($6, cigarSegs, /[0-9]+[MIDNSHPX=]/)
            endOfSpliter2InR2 = 0
            for (i = 1; i <= n; ++i)
            {
                patsplit(cigarSegs[i], num, /[0-9]+/, labels)
                if (labels[1] ~ /[MI=X]/ || labels[1] ~ /[SH]/ && i == 1)
                    endOfSpliter2InR2 += num[1]
            }
            print $2, $3, endOfSpliter2InR2
        }
    }
    '
}

appendSpliter2SeqRef12()
{
    # Usage: appendSpliter2SeqRef12 <R2map
    # Append the sequence of spliter2, as well as ref1 and ref2 and sgRNA
    # Input: flag|spliter2|alignmentEndInQuery
    # Output: flag|spliter2|alignmentEndInQuery|spilter2seq|ref1|ref2|sgRNA
    awk -F "\t" -v spliter2fa="${spliter2}.fa" -v ref12="$ref12" -v sgRNAfile=${sgRNAfile} -v OFS="\t" '
    BEGIN{
        while (getline spliter2name <spliter2fa)
        {
            getline spliter2seq <spliter2fa
            name = substr(spliter2name, 2)
            name2seq[name] = spliter2seq
            getline r12 <ref12
            name2r12[name] = r12
            getline sgRNA <sgRNAfile
            name2sgRNA[name] = sgRNA
        }
    }
    {
        print $1, $2, $3, toupper(name2seq[$2]), toupper(name2r12[$2]), toupper(name2sgRNA[$2])
    }
    '
}

R1map()
{
    # Usage: R1map <R1reads
    # Map R1 to spliter1, and retrieve flag and matched spliter1
    # Input: R1reads
    # Output: flag|spilter1
    bowtie2 --quiet --norc --mm --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,$minscoreR1 -r -x "$spliter1" -U - 2>/dev/null | samtools view | cut -f2,3
}

FilterSpilters()
{
    # Usage: FilterSpilters <R2|R1|count|R2CutAdapt|flag2|spliter2|endOfSpliter2InR2|spliter2seq|ref1|ref2|sgRNA|flag1|spliter1
    # Input: R2|R1|count|R2CutAdapt|flag2|spliter2|endOfSpliter2InR2|spliter2seq|ref1|ref2|sgRNA|flag1|spliter1
    # Output: R2|count|R2CutAdapt|endOfSpliter2InR2|spliter2seq|ref1|ref2|sgRNA
    # $1:R2|$3:count|$4:R2CutAdapt|$5:flag2|$6:spliter2|$7:endOfSpliter2InR2|$8:spliter2seq|$9:ref1|$10:ref2|$11:sgRNA|$12:flag1|$13:spliter1
    # filter the read pair if one of the following happens:
    # 1. R2 does not match any record in spliter2
    # 2. R1 does not match any record in spliter1
    # 3. spliter1 matching R1 and spliter2 matching R2 are not paired
    # 4. After removing 5' spliter2 and 3' adapter, R2 is shorter than minQueryR2
    awk -F "\t" -v OFS="\t" -v minQueryR2=$minQueryR2 -v fqR1="$fqR1" '
    {
        if (($5/4)%2 == 1 || ($12/4)%2 == 1 || $6 != $13 || $7 + minQueryR2 > length($4))
            print $0 > fqR1 ".not_find";
        else
        {
            print $1, $3, $4, $7, $8, $9, $10, $11;
            totalCount += $3
        }
    }
    END{
        print totalCount > fqR1 ".total";
    }
    ' 
}

cumulate_R2_sort()
{
    # Usage: cumulate_R2 <R2|count|R2CutAdapt|endOfSpliter2InR2|spliter2seq|ref1|ref2|sgRNA
    # Input: R2|count|R2CutAdapt|endOfSpliter2InR2|spliter2seq|ref1|ref2|sgRNA
    # Output: R2|count|R2CutAdapt|endOfSpliter2InR2|spliter2seq|ref1|ref2|sgRNA
    # Cumulate the adjacent duplicate R2 count. Sort the result first by the dict order of spliter2seq, and then by the descent order of count
    awk -F "\t" -v OFS="\t" '
    {
        if ($1 != R2)
        {
            if (NR > 1)
                print R2, count, R2CutAdapt, endOfSpliter2InR2, spliter2seq, ref1, ref2, sgRNA
            R2 = $1
            count = $2
            R2CutAdapt = $3
            endOfSpliter2InR2 = $4
            spliter2seq = $5
            ref1 = $6
            ref2 = $7
            sgRNA = $8
        }
        else
            count += $2
    }
    END{
        print R2, count, R2CutAdapt, endOfSpliter2InR2, spliter2seq, ref1, ref2, sgRNA
    }
    ' | sort -k5,5 -k2,2nr
}

fqR1=$1
fqR2=$2
spliter1=$3
spliter2=$4
sgRNAfile=$5
ref12=$6
minscoreR1=$7
minscoreR2=$8
minQueryR2=$9
project_path="$(dirname $(realpath $0))/../.."

# list R2 and R1 side-by-side, and count duplicates: R2\tR1\tcount
$project_path/pv-1.8.5/pv -c -N "count $fqR1" "$fqR2" | sed -n '2~4p' | paste - <(sed -n '2~4p' "$fqR1") | sort | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' >"$fqR1.count"

$project_path/pv-1.8.5/pv -c -N "demultiplex $fqR1" "$fqR1.count" | cut -f1 | R2map | appendSpliter2SeqRef12 | paste "$fqR1.count" <(cut -f1 "$fqR1.count" | cutadaptPlain GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC) - <(cut -f2 "$fqR1.count" | R1map) | FilterSpilters | cumulate_R2_sort >"${fqR1}.demultiplex"