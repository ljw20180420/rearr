#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage

```shell
$ getSxPlasmidFileRef.sh \
    plasmid_file \
    genome \
    bowtie2index \
    [ext1up ext1down ext2up ext2down]
```

<pre class="mermaid">
---
title: getSxPlasmidFileRef
---
flowchart TD
    PF[(
        <h1>plasmid_file</h1>
        <table>
            <tr>
                <th>adapter#0040;20bp#0041; + sgRNA#0040;20bp#0041; + scaffold#0040;83/93bp#0041; + query#0040;44bp#0041; + 3bp + RCbarcode#0040;18bp#0041; + RCprimer#0040;21bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )] --> GSPFR[getSxPlasmidFileRef.sh]
    GF[(<h1>genome_file</h1>)] --> GSPFR
    BI[(<h1>bowtie2index</h1>)] --> GSPFR
    EXT[(
        <h1>extentions</h1>
        <table>
            <tr>
                <td>ext1up</td>
                <td>ext1down</td>
                <td>ext2up</td>
                <td>ext2down</td>
            </tr>
        </table>
    )] --> GSPFR
    GSPFR --> REF[(
        <h1>stdout</h1>
        <table>
            <tr>
                <th>start1</th>
                <th>ref1</th>
                <th>end1</th>
                <th>start2</th>
                <th>ref2</th>
                <th>end2</th>
            </tr>
            <tr>
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

- Extract references from the in-house `plasmid_file` and output to `stdout`.
- For `NGG` `plasmid_file`, the 44bp `query` can be perfectly mapped to the genome.
- For `NAA` `plasmid_file`, 17~18bp of `query` is `TT`, which should be replaced by `CC` in order to map genome. 
- The actual cut site is inferred from mapping result.
- `ref1` consists of `ext1up` bases upstream to the cut site and `ext1down` bases downstream to the cut site.
- `ref2` consists of `ext2up` bases upstream to the cut site and `ext2down` bases downstream to the cut site.
- Note that for `NAA` `plasmid_file`, the retrieved reference need replace `GG` (`query` and reference have opposite strands, so `CC` becomes `GG`) back to `AA`.

# Source
~~~bash
plasmid_file=$1
genome=$2
bowtie2index=$3
ext1up=${4:-50}
ext1down=${5:-0}
ext2up=${6:-10}
ext2down=${7:-100}


getSxPlasmidFileTarget.pl "${plasmid_file}" | bowtie2 --quiet --mm -x "${bowtie2index}" -r -U - 2> /dev/null | samtools view | gawk -f sxTargetSam2Bed.awk -- ${ext1up} ${ext1down} ${ext2up} ${ext2down} | bedtools getfasta -s -fi "${genome}" -bed - | sed '1~2d' | getSxRefFile.pl ${ext1up} ${ext2up} ${plasmid_file: -6:1}
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
