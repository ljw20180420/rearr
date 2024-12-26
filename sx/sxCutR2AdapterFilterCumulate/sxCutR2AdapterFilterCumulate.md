#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
sxCutR2AdapterFilterCumulate.md demultiplexFile minToMapShear >toMapFile
```

# Introduction
The output of [`demultiplex.md`][`demultiplex.md`] needs further post-process before feed to [`rearr`][`rearr`]. For Shi Xing's data, this is done by this in-house script. The output format fits the input format of [`rearr`][`rearr`].
```bash
query<tab>count<tab>refId<newline>
```

The post-process consists of three steps.
1. Remove `adapter` from 3' of R2.
2. Remove `primer`, `barcode` and a 3bp gap from 5' of R2.
3. Filter out if the remain of R2 is shorter than minToMapShear.
4. Accumulate the adjacent duplicates by [`sxCumulateToMapCutAdaptSpliter.awk`][`sxCumulateToMapCutAdaptSpliter.awk`].

[`demultiplex.md`]: /sx_lcy/core/demultiplex/
[`rearr`]: /sx_lcy/core/rearr/
[`sxCumulateToMapCutAdaptSpliter.awk`]: /sx_lcy/sx/sx-cut-r2-adapter-filter-cumulate/sx-cumulate-to-map-cut-adapt-spliter/

# Source
~~~bash
cutadaptPlain()
{
    # Usage: cutadaptPlain <plainseq 3'adapter
    # cutadapt does not accept plainseq. This function transform plainseq to fasta before feed to cutadapt, and then transform the fasta output back to plainseq
    # Input: plainseq
    # Output: 3' trimmed plainseq
    sed '=' | sed '1~2s/^/>s/' | cutadapt -a $1 - 2> /dev/null | sed '1~2d'
}

rmDupFile=$1
minToMapShear=$2
cut -f1 $rmDupFile | cutadaptPlain GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC | paste - <(cut -f3-4,8 $rmDupFile) | gawk -F "\t" -v OFS="\t" -v minToMapShear=$minToMapShear '
{
    if ($4 + minToMapShear <= length($1)) {
        print substr($1, $4 + 1), $2, $3
    }
}' | gawk -f sxCumulateToMapCutAdaptSpliter.awk
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
