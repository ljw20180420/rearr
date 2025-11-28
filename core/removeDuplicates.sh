#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage

```bash
$ removeDuplicates.sh fastqR1 [fastqR2 [fastqR3 ...]]
```

<pre class="mermaid">
---
title: removeDuplicates.sh
---
flowchart TD
    R1[(fastqR1)] --> RD[removeDuplicates.sh]
    R2[(faqstR2)] --> RD
    RN[(...)] --> RD
    RD --> UNIQUE[(
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
</pre>

- The input fastq files can be gz or zip compressed.
- Multiple fastq files, even more than two, are supported in theory. However, in practice, one has only `R1` and `R2`.
- Two sequencings are duplicates if they are the same across all fastq file. The duplication number `#` follows the tandom tab separated list of reads in fastq files in `stdout` of `removeDuplicates.sh`.

# Why use several `fastq` files as input of `removeDuplicates.sh`
The paired-end next-generation sequencing (NGS) is quite common. Although mappable segment may be only in `R1` or `R2`, the other end still helps to determine the locus of the sequence. See [`demultiplex.sh`][demultiplex.sh.md].

# Should I directly input raw `fastq` file, or remove `adapter`, `barcode` and so on before the input into `removeDuplicates.sh`
The `stdout` of `removeDuplicates.sh` are aligned to the so-call `spliters` in [`demultiplex.sh`][demultiplex.sh.md] to determine the loci of lines. If you preserve `adapter`, `barcode` and so on in the input `fastq` files, it is suggested to provide them in `spliters` as well.

[demultiplex.sh.md]: /auxilary/demultiplex.sh.html

# Source
~~~bash
fqlist=""
for fq in "$@"
do
    if (file $fq | grep -q compressed)
    then
        fqlist="$fqlist <(zcat $fq)"
    else
        fqlist="$fqlist $fq"
    fi
done

eval paste $fqlist | sed -n '2~4p' | sort | uniq -c | gawk '
    {
        for (i = 2; i <= NF; ++i)
            printf("%s\t", $i)
        print $1
    }
'
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
