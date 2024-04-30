#!/bin/bash
# Usage: demultiplex.sh inputFile spliterToMap spliterPair minScoreToMap minScorePair
# spliter(ToMap|Pair) is a bowtie2 index used to split (toMap|pair)
# minScore(ToMap|Pair) thres the match between (toMap|pair) and spliter(ToMap|Pair)

# inputFile format
# toMap|pair|count
# toMap is the read to map (either R1 or R2)
# pair is the read paired with toMap
# count is the duplicate number of (toMap, pair)

mapSpliter()
{
    # Usage: mapSpliter minScore spliterBowtie2Index <reads
    # Map reads to spliter by bowtie2
    # Input: plain reads
    # Output: sam file without header
    minScore=$1
    spliter=$2
    bowtie2 --quiet --norc --mm --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,${minScore} -r -x "${spliter}" -U - 2>/dev/null | samtools view
}

FilterSpilters()
{
    # Input: toMap|Pair|count|flagToMap|spliterToMap|endOfSpliterPosToMap|flagPair|spliterPair|endOfSpliterPosPair
    # Output: toMap|Pair|count|spliter|endOfSpliterPosToMap|endOfSpliterPosPair
    # filter row if one of the following happens:
    # 1. toMap does not match any record in spliterToMap
    # 2. pair does not match any record in spliterPair
    # 3. spliterToMap and spliterPair do not match
    gawk -F "\t" -v OFS="\t" '
    {
        if (($4/4)%2 == 0 && ($7/4)%2 == 0 && $5 == $8)
            print $1, $2, $3, $5, $6, $9
    }
    ' 
}

inputFile=$1
spliterToMap=$2
spliterPair=$3
minScoreToMap=$4
minScorePair=$5


pv -c -N "demultiplex ${inputFile}" "${inputFile}" | cut -f1 | mapSpliter ${minScoreToMap} ${spliterToMap} | gawk -f endOfSpliterPos.awk | paste "${inputFile}" - <(cut -f2 "${inputFile}" | mapSpliter ${minScorePair} ${spliterPair} | gawk -f endOfSpliterPos.awk) | FilterSpilters >"${inputFile}.demultiplex"