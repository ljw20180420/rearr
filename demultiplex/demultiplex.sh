#!/bin/bash
# Usage: spliterTarget=spliterTarget spliterPair=spliterPair minScoreTarget=minScoreTarget minScorePair=minScorePair demultiplex.sh inputFile >demultiplexFile
# spliter(Target|Pair) is a bowtie2 index used to split (Target|pair)
# minScore(Target|Pair) thres the match between (Target|pair) and spliter(Target|Pair)

# inputFile format
# Target|pair|count
# Target is the read to map (either R1 or R2)
# pair is the read paired with Target
# count is the duplicate number of (Target, pair)

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
    if [ -n "$spliterTarget" ] && [ -n "$spliterPair" ]
    then
        # both spliterTarget and spliterPair are specified
        # Input: Target|Pair|count|flagTarget|spliterTarget|endOfSpliterPosTarget|flagPair|spliterPair|endOfSpliterPosPair
        # Output: Target|Pair|count|spliter|endOfSpliterPosTarget|endOfSpliterPosPair
        # filter row if one of the following happens:
        # 1. Target does not match any record in spliterTarget
        # 2. pair does not match any record in spliterPair
        # 3. spliterTarget and spliterPair do not match
        gawk -F "\t" -v OFS="\t" '
        {
            if (($4/4)%2 == 0 && ($7/4)%2 == 0 && $5 == $8)
                print $1, $2, $3, $5, $6, $9
        }
        '
    else
        # only spliter(Target|Pair) is specified
        # Input: Target|Pair|count|flag(Target|Pair)|spliter(Target|Pair)|endOfSpliterPos(Target|Pair)
        # Output: Target|Pair|count|spliter|endOfSpliterPos(Target|Pair)
        # filter row if one of the following happens:
        # 1. (Target|Pair) does not match any record in spliter(Target|Pair)
        gawk -F "\t" -v OFS="\t" '
        {
            if (($4/4)%2 == 0)
                print $1, $2, $3, $5, $6
        }
        '
    fi
}

inputFile=$1
if [ -z "$spliterTarget" ] && [ -z "$spliterPair" ]
then
    echo "At least one of spliterTarget or spliterPair need to be specified."
fi
if [ -n "$spliterTarget" ] && [ -z "$minScoreTarget" ]
then
    echo "spliterTarget is specified but minScoreTarget is not."
fi
if [ -n "$spliterPair" ] && [ -z "$minScorePair" ]
then
    echo "spliterPair is specified but minScorePair is not."
fi

map1=""
if [ -n "$spliterTarget" ] 
then
    map1="<(cut -f1 ${inputFile} | mapSpliter ${minScoreTarget} ${spliterTarget} | gawk -f endOfSpliterPos.awk)"
fi
map2=""
if [ -n "$spliterPair" ] 
then
    map2="<(cut -f2 "${inputFile}" | mapSpliter ${minScorePair} ${spliterPair} | gawk -f endOfSpliterPos.awk)"
fi

eval paste "${inputFile}" $map1 $map2 | FilterSpilters