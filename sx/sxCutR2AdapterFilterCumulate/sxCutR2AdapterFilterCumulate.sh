#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
sxCutR2AdapterFilterCumulate.sh demultiplex_file minToMapShear
```

<pre class="mermaid">
---
title: sxCutR2AdapterFilterCumulate.sh
---
flowchart TD
    ONTARGET[(
        <h1>demultiplex_file</h1>
        <table>
            <tr>
                <th>R2</th>
                <th>R1</th>
                <th>#</th>
                <th>id</th>
                <th>rstart2</th>
                <th>rend2</th>
                <th>qstart2</th>
                <th>qend2</th>
                <th>rstart1</th>
                <th>rend1</th>
                <th>qstart1</th>
                <th>qend1</th>
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
    )] --> sxCRAFC[sxCutR2AdapterFilterCumulate.sh]
    sxCRAFC --> QUERY[(
        <h1>stdout</h1>
        <table>
            <tr>
                <th>query</th>
                <th>#</th>
                <th>id</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]
</pre>

- The `stdout` of [`demultiplex.sh`][demultiplex.sh.md] (`demultiplex_file`) needs further post-process before feed to [core part of rearr][core.md]. Note that we put `R2` before `R1` for in-house data when calling [`demultiplex.sh`][demultiplex.sh.md].
- Our in-house data have `query` on R2.
  > `R2` = `primer`(21bp) + `barcode`(18bp) + 3bp + `query`(44bp) + `RCscaffold`(83/93bp).
  We extract `query` from `R2` as follows.
    1. Trim 3' `RCscaffold`.
    2. Extract `query` 3bp downstream to `qend2` (the alignment end position of `barcode`).
    3. Filter out `query` shorter than `minToMapShear`.
    4. Accumulate the adjacent duplicates of `query` by [`sxCumulateToMapCutAdaptSpliter.awk`][sxCumulateToMapCutAdaptSpliter.awk.md].  

[demultiplex.sh.md]: /auxilary/demultiplex.sh.html
[core.md]: /core.html
[sxCumulateToMapCutAdaptSpliter.awk.md]: /sx/sxCutR2AdapterFilterCumulate/sxCumulateToMapCutAdaptSpliter.awk.html

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
    if ($4 + 3 + minToMapShear <= length($1)) {
        print substr($1, $4 + 4), $2, $3
    }
}' | gawk -f sxCumulateToMapCutAdaptSpliter.awk
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
