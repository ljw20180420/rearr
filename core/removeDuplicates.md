#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
removeDuplicates.md fastq1(.gz|.zip) fastq2(.gz|.zip) fastq3(.gz|.zip) ... >rmDupFile
```

# Introduction
The input of `removeDuplicates.md` are several `fastq` files. `removeDuplicates.md` accepts `fastq` files in `.gz` or `.zip` compressed form. The `Nth` sequence of each input `fastq` file forms the `Nth` record. Two records are deplicates if their component sequences from the same input `fastq` file are always the same. `removeDuplicates.md` remove and count duplicated records. The `stdout` are lines of the form
```
seq1<tab>seq2<tab>...<tab>count<newline>
```

## Why use several `fastq` files as input of `removeDuplicates.md`
The paired-end next-generation sequencing (NGS) is quite common. Although mappable segment may be only in `R1` or `R2`, the other end still helps to determine the locus of the sequence. See [`demultiplex.md`][`demultiplex.md`].

## Should I directly input raw `fastq` file, or remove `adapter`, `barcode` and so on before the input into `removeDuplicates.md`
The `stdout` of `removeDuplicates.md` are aligned to the so-call `spliters` in [`demultiplex.md`][`demultiplex.md`] to determine the loci of lines. If you preserve `adapter`, `barcode` and so on in the input `fastq` files, it is suggested to provide them in `spliters` as well.

[`demultiplex.md`]: /rearr/core/demultiplex/

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
