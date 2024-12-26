#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
getSxCsvFileRef.md csvfile genome bowtie2index [ext1up ext1down ext2up ext2down]
```

# Introcution
This is an in-house script to extract references from the `csvfile` of sx and lcy. The composition of the input csvfile is 
```
adapter(20bp) + sgRNA(20bp) + scaffold(83/93bp) + target(44bp) + 3bp + RCbarcode(18bp) + RCprimer(21bp)
```
For `NGG` csvfile, the 44bp `target` can be perfectly mapped to the genome. For `NAA` csvfile, 17~18bp of target is `TT`, which should be replaced by `CC` in order to map genome. After mapping, the actual cut site is inferred. `ref1` consists of `ext1up` bases upstream to the cut site and `ext1down` bases downstream to the cut site. `ref2` is composed similarly. Note that for `NAA` csvfile, the retrieved reference need replace `GG` (target and reference always have opposite strands, so `CC` becomes `GG`) back to `AA`.

# Source
~~~bash
csvfile=$1
genome=$2
bowtie2index=$3
ext1up=${4:-50}
ext1down=${5:-0}
ext2up=${6:-10}
ext2down=${7:-100}


getSxCsvFileTarget.pl "${csvfile}" | bowtie2 --quiet --mm -x "${bowtie2index}" -r -U - 2> /dev/null | samtools view | gawk -f sxTargetSam2Bed.awk -- ${ext1up} ${ext1down} ${ext2up} ${ext2down} | bedtools getfasta -s -fi "${genome}" -bed - | sed '1~2d' | getSxRefFile.pl ${ext1up} ${ext2up} ${csvfile: -6:1}
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
