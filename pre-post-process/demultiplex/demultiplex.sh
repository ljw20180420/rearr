#!/bin/bash
# Usage: spliterIndices=index1,index2,... minScores=score1,score2,... demultiplex.sh inputFile >demultiplexFile
# spliterIndices are bowtie2 index aligned by sequences in inputFile (i.e. the stdout of removeDuplicates.sh). minScores are feed to bowtie2 to filter alignments with low scores. Note that the order of spliterIndices and minScores matters. The alignments are in local mode instead of end-to-end mode, and there is no reverse complement. The full alignment setting is --norc --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,scoreN.

mapSpliter()
{
    # Usage: mapSpliter minScore spliterIndex <reads
    # Map reads to spliter by bowtie2
    # Input: plain reads
    # Output: sam file without header
    minScore=$1
    spliterIndex=$2
    bowtie2 --quiet --mm --norc --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,${minScore} -r -x "${spliterIndex}" -U - 2>/dev/null | samtools view
}

filterSpilters()
{
    # Input: seq1|seq2|...|count|flag1|spliter1|spliterStart1|spliterEnd1|seqStart1|seqEnd1|flag2|spliter2|spliterStart2|spliterEnd2|seqStart2|seqEnd2|...
    # Output: seq1|seq2|...|count|spliter|spliterStart1|spliterEnd1|seqStart1|seqEnd1|spliterStart2|spliterEnd2|seqStart2|seqEnd2
    # filter out rows with one of the following happens:
    # 1. seqN failed to align spliterIndexN
    # 2. spilterM != spliterN for some M and N
    gawk -F "\t" -v OFS="\t" -v firstFlagPos=$((${#spliterIndexArray[@]}+2)) '
        {
            spliter = $(firstFlagPos+1)
            for (i = firstFlagPos; i <= NF; i += 6) {
                if ($i != 0 || $(i+1) != spliter) {
                    next
                }
            }
            for (i = 1; i < firstFlagPos; ++i) {
                printf("%s\t", $i)
            }
            printf("%s\t", $(firstFlagPos+1))
            for (i = firstFlagPos; i < NF; ++i) {
                if ((i - firstFlagPos) % 6 <= 1) {
                    continue
                }
                printf("%s\t", $i)
            }
            printf("%s\n", $NF)
        }
    '
}

inputFile=$1

IFS=',' read -r -a spliterIndexArray <<< "$spliterIndices"
IFS=',' read -r -a minScoreArray <<< "$minScores"
maps=""
for ii in ${!spliterIndexArray[@]}
do
    maps="${maps} <(cut -f$((ii+1)) ${inputFile} | mapSpliter ${minScoreArray[$ii]} ${spliterIndexArray[$ii]} | gawk -f getAlignPos.awk)"
done

eval paste "${inputFile}" $maps | filterSpilters