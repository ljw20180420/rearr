#!/bin/bash
# Usage: demultiplex.sh inputFile spliterToMap spliterPair minScoreToMap minScorePair minToMapShear
# spliter(ToMap|Pair) is a bowtie2 index used to split (toMap|pair)
# minScore(ToMap|Pair) thres the match between (toMap|pair) and spliter(ToMap|Pair)
# Remove spliterToMap from 5' and adapter from 3' of toMap. The remain must contain at least minToMapShear bases.

# inputFile format
# toMap|pair|count
# toMap is the read to map (either R1 or R2)
# pair is the read paired with toMap
# count is the duplicate number of (toMap, pair)

cutadaptPlain()
{
    # Usage: cutadaptPlain <plainseq 3'adapter
    # cutadapt does not accept plainseq. This function transform plainseq to fasta before feed to cutadapt, and then transform the fasta output back to plainseq
    # Input: plainseq
    # Output: 3' trimmed plainseq
    sed '=' | sed '1~2s/^/>s/' | cutadapt -a $1 - 2> /dev/null | sed '1~2d'
}

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
    # Usage: FilterSpilters <toMapCutAdapt|count|flagToMap|spliterToMap|endOfSpliterPos|flagPair|spliterPair
    # Input: toMapCutAdapt|count|flagToMap|spliterToMap|endOfSpliterPos|flagPair|spliterPair
    # Output: toMapCutAdaptSpliter|count|spliter
    # filter row if one of the following happens:
    # 1. toMap does not match any record in spliterToMap
    # 2. pair does not match any record in spliterPair
    # 3. spliterToMap and spliterPair do not match
    # 4. After removing 5' spliterToMap from toMapCutAdapt, the resulting toMapCutAdaptSpliter is shorter than minToMapShear
    gawk -F "\t" -v OFS="\t" -v minToMapShear=$1 '
    {
        if (($3/4)%2 == 0 && ($6/4)%2 == 0 && $4 == $7 && $5 + minToMapShear <= length($1))
            print substr($1, $5 + 1), $2, $4
    }
    ' 
}

inputFile=$1
spliterToMap=$2
spliterPair=$3
minScoreToMap=$4
minScorePair=$5
minToMapShear=$6

pv -c -N "demultiplex ${inputFile}" "${inputFile}" | cut -f1 | mapSpliter ${minScoreToMap} ${spliterToMap} | gawk -f endOfSpliterPos.awk | paste <(cut -f1 "${inputFile}" | cutadaptPlain GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC) <(cut -f3 "${inputFile}") - <(cut -f2 "${inputFile}" | mapSpliter ${minScorePair} ${spliterPair} | cut -f2,3) | FilterSpilters ${minToMapShear} | gawk -f cumulateToMapCutAdaptSpliter.awk >"${inputFile}.demultiplex"