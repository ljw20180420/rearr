#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage

```bash
$ markerIndices=marker1,marker2,... \
  minScores=score1,score2,... \
  demultiplex.sh \
    removeDuplicates_file
```

<pre class="mermaid">
---
title: demultiplex.sh
---
flowchart TD
    UNIQUE[(
        <h1>removeDuplicates_file</h1>
        <table>
            <tr>
                <th>R1</th>
                <th>R2</th>
                <th>...</th>
                <th>#</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]
    UNIQUE --> DM[demultiplex.sh]
    MK1[(marker1)] --> DM
    MK2[(marker2)] --> DM
    MKN[(...)] --> DM
    SCORE[(
        <h1>minScores</h1>
        <table>
            <tr>
                <th>score1</th>
                <th>score2</th>
                <th>...</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]
    SCORE --> DM
    DM --> ONTARGET[(
        <h1>stdout</h1>
        <table>
            <tr>
                <th>R1</th>
                <th>R2</th>
                <th>...</th>
                <th>#</th>
                <th>id</th>
                <th>rstart1</th>
                <th>rend1</th>
                <th>qstart1</th>
                <th>qend1</th>
                <th>rstart2</th>
                <th>rend2</th>
                <th>qstart2</th>
                <th>qend2</th>
                <td>...</td>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]
</pre>

- `markerN` is `bowtie2` index of a `fasta` file containing possible alignment targets of `RN`.
- `scoreN` is feed to `bowtie2` to filter alignments with low scores. The alignments are in local mode instead of end-to-end mode, and there is no reverse complement. The full alignment setting is
```bash
--norc --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,scoreN
```
- The fasta name of `marker1` and `marker2` must the 0-based reference `id`, see the [core part of rearr][core.md].
- `R1` and `R2` are given in `removeDuplicates_file`, which is `stdout` of [`removeDuplicates.sh`][removeDuplicates.sh.md].
- `demultiplex.sh` only output a row in `removeDuplicates_file` only of both `R1` and `R2` aligned to targets with consistant reference `id`.
- `R1`, `R2` and `#` are copied from `removeDuplicates_file`.
- `rstartN` and `rendN` denotes the left-close-right-open 0-based range of the aligned part of the target of `RN`. `qstartN` and `qendN` denotes that of `RN`.
- Multiple `marker` (more than two) are supported in theory.

[core.md]: /core.html
[removeDuplicates.sh.md]: /auxilary/removeDuplicates.sh.html

# Source
~~~bash
mapMarker()
{
    # Usage: mapMarker minScore markerIndex <reads
    # Map reads to marker by bowtie2
    # Input: plain reads
    # Output: sam file without header
    minScore=$1
    markerIndex=$2
    bowtie2 --quiet --mm --norc --local -L 15 --ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1 --score-min C,${minScore} -r -x "${markerIndex}" -U - 2>/dev/null | samtools view
}

filterMarkers()
{
    # Input: seq1|seq2|...|count|flag1|refId1|markerStart1|markerEnd1|seqStart1|seqEnd1|flag2|refId2|markerStart2|markerEnd2|seqStart2|seqEnd2|...
    # Output: seq1|seq2|...|count|refId|markerStart1|markerEnd1|seqStart1|seqEnd1|markerStart2|markerEnd2|seqStart2|seqEnd2|...
    # filter out rows with one of the following happens:
    # 1. seqN failed to align markerIndexN
    # 2. refIdM != refIdN for some M and N
    gawk -F "\t" -v OFS="\t" -v firstFlagPos=$((${#markerIndexArray[@]}+2)) '
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

IFS=',' read -r -a markerIndexArray <<< "${markerIndices}"
IFS=',' read -r -a minScoreArray <<< "$minScores"
maps=""
for ii in ${!markerIndexArray[@]}
do
    maps="${maps} <(cut -f$((ii+1)) ${inputFile} | mapMarker ${minScoreArray[$ii]} ${markerIndexArray[$ii]} | gawk -f getAlignPos.awk)"
done

eval paste "${inputFile}" $maps | filterMarkers
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
