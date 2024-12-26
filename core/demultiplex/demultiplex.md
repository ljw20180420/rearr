#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
spliterIndices=index1,index2,... minScores=score1,score2,... demultiplex.md inputFile >demultiplexFile
```

# Introduction
Each `spliterIndice` is `bowtie2` index of a `fasta` file containing possible alignment references (i.e. `spliter`s) of the corresponding (depends on the order) output `seq` of [`removeDuplicates.md`][`removeDuplicates.md`]. The `Nth` `spliter`s of `spliterIndice`s are usually different. However, they must have the same name. The name denotes the 0-based reference id (`refId`).

`bowtie2` either aligns `seq` to one of the `spliter`s or failed to align. If all `seq`s in a line output of [`removeDuplicates.md`][`removeDuplicates.md`] successfully align to `spliter`s of the same `refId`, then `demultiplex.md` will print the following to `stdout`.
```
seq1<tab>seq2<tab>...<tab>count<tab>refId<tab>rs1<tab>re1<tab>qs1|qe1|rs2|re2|qs2|qe2|...
```
`seqN` and `count` are copied from `inputFile`. `rsN` and `reN` denotes the left-close-right-open 0-based range of the aligned part of the `Nth` reference (`spliter`) in the local alignment. `qsN` and `qeN` denotes that of the query (`seq`).

## Alignment details
`minScore`s are feed to `bowtie2` to filter alignments with low scores. The alignments are in local mode instead of end-to-end mode, and there is no reverse complement. The full alignment setting is
```bash
--norc --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,scoreN
```

[`removeDuplicates.md`]: /sx_lcy/core/remove-duplicates/

# Source
~~~bash
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
    # Input: seq1|seq2|...|count|flag1|refId1|spliterStart1|spliterEnd1|seqStart1|seqEnd1|flag2|refId2|spliterStart2|spliterEnd2|seqStart2|seqEnd2|...
    # Output: seq1|seq2|...|count|refId|spliterStart1|spliterEnd1|seqStart1|seqEnd1|spliterStart2|spliterEnd2|seqStart2|seqEnd2|...
    # filter out rows with one of the following happens:
    # 1. seqN failed to align spliterIndexN
    # 2. refIdM != refIdN for some M and N
    gawk -F "\t" -v OFS="\t" -v firstFlagPos=$((${#spliterIndexArray[@]}+2)) '
        {
            refId = $(firstFlagPos+1)
            for (i = firstFlagPos; i <= NF; i += 6) {
                if ($i != 0 || $(i+1) != refId) {
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
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
