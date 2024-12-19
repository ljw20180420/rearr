# `removeDuplicates.sh`
The input of `removeDuplicates.sh` are several `fastq` files. `removeDuplicates.sh` accepts `fastq` files in `.gz` or `.zip` compressed form. The `Nth` sequence of each input `fastq` file forms the `Nth` record. Two records are deplicates if their component sequences from the same input `fastq` file are always the same. `removeDuplicates.sh` remove and count duplicated records. The `stdout` are lines of the form
```
seq1<tab>seq2<tab>...<tab>count<newline>
```

## Why use several `fastq` files as input of `removeDuplicates.sh`
The paired-end next-generation sequencing (NGS) is quite common. Although mappable segment may be only in `R1` or `R2`, the other end still helps to determine the locus of the sequence. See [`demultiplex.sh`](#demultiplex).

## Should I directly input raw `fastq` file, or remove `adapter`, `barcode` and so on before the input into `removeDuplicates.sh`
The `stdout` of `removeDuplicates.sh` are aligned to the so-call `spliters` in [`demultiplex.sh`](#demultiplex) to determine the loci of lines. If you preserve `adapter`, `barcode` and so on in the input `fastq` files, it is suggested to provide them in `spliters` as well.

## Usage
```bash
removeDuplicates.sh fastq1(.gz|.zip) fastq2(.gz|.zip) fastq3(.gz|.zip) ... >rmDupFile
```


<a id="demultiplex"></a>
# `demultiplex.sh`


# Usage: spliterTarget=spliterTarget spliterPair=spliterPair minScoreTarget=minScoreTarget minScorePair=minScorePair demultiplex.sh inputFile >demultiplexFile
# spliter(Target|Pair) is a bowtie2 index used to split (Target|pair)
# minScore(Target|Pair) thres the match between (Target|pair) and spliter(Target|Pair)

# inputFile format
# Target|pair|count
# Target is the read to map (either R1 or R2)
# pair is the read paired with Target
# count is the duplicate number of (Target, pair)

The `Nth` sequence in an `stdout` line, say `seqN`, is