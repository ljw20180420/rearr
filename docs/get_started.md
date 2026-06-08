---
title: Get started
---

# Example

Here is a minimal example to run [`rearr`][rearr] on sgRNA library.

- [`example.fq`][example.fq] contains the sequencing reads.
- [`example.fa`][example.fa] contains markers. These markers appear in sequences in [`example.fq`][example.fq]. They are used to demultiplex sequencing reads to the correct reference.
- [`example.ref`][example.ref] contains references. Each reference (row) in [`example.ref`][example.ref] corresponds a marker in [`example.fa`][example.fa]. If a read hits, say the 5th marker in [`example.fa`][example.fa], then this reads will be mapped to the reference in the 5th row in [`example.ref`][example.ref].

## Structure of references in [`example.ref`][example.ref]

Each row in [`example.ref`][example.ref] has the following form.

```
┌───┬────────────────────┬──────────────────────────────────────────────────┬──────────────────────────────────────────────────────┬──────────────────────┬─────────────────────────────┐
│ 0 | upstream reference | upstream cleavage position in upstream reference | downstream cleavage position in downstream reference | downstream reference | downstream reference length |
└───┴────────────────────┴──────────────────────────────────────────────────┴──────────────────────────────────────────────────────┴──────────────────────┴─────────────────────────────┘
```

[`rearr`][rearr] use this format to handle double cleavage CRISPR/Cas protocols. Note that single cleavage protocols are just special cases of double cleavage protocols. Suppose you design a pair of sgRNA to cut two sites on the genome.

```
site1:TTTTTCACATTAAAGTAGCTGGATTGTACCTAAATTCATCCGTACGGGAT|ACTGGTGCGGTTCGTATCGGGAACGAGCTACACCGCTCTATCGGGCCCGG
site2:TCACTTTCTAACGCCCCACCCGACTGTATTAGATTTGCCGTTCAAGAGGC|AAAAATGAAGGGAGATAGTTTTCAGCAGGAACGGAGCTGAGGTAAGGCTA
```

You want the oberve the indels of the junction between the upstream of `site1` and downstream of `site2`.

```
TTTTTCACATTAAAGTAGCTGGATTGTACCTAAATTCATCCGTACGGGAT|AAAAATGAAGGGAGATAGTTTTCAGCAGGAACGGAGCTGAGGTAAGGCTA
```

In this case,

```
~~~~~~~~~~~~~~~~~~~~~TTTTTCACATTAAAGTAGCTGGATTGTACCTAAATTCATCCGTACGGGAT|ACTGGTGCGGTTCGTATCGGGAACGAGCTACACCGCTCTATCGGGCCCGG
  upstream reference:TTTTTCACATTAAAGTAGCTGGATTGTACCTAAATTCATCCGTACGGGAT|ACTGGTGCGG
                     TCACTTTCTAACGCCCCACCCGACTGTATTAGATTTGCCGTTCAAGAGGC|AAAAATGAAGGGAGATAGTTTTCAGCAGGAACGGAGCTGAGGTAAGGCTA
downstream reference:                                        TTCAAGAGGC|AAAAATGAAGGGAGATAGTTTTCAGCAGGAACGGAGCTGAGGTAAGGCTA
```

`upstream reference` covers a large part upstream to the cleavage site as well as a small part downstream to the cleavage site. The small part downstream to the cleavage site is used to capture predictable insertions. Similarly, `downstream reference` covers a large part downstream to the cleavage site as well as a small part upstream to the cleavage site. The small part upstream to the cleavage site is used to capture predictable insertions. `upstream cleavage position in upstream reference` is the length of the part of `upstream reference` upstream to the cleavage site (50 in this example). `downstream cleavage position in downstream reference` is the length of the part of `downstream reference` upstream to the cleavage site (10 in this example). You may adjust the size of small part across the cleavage site in `upstream reference` and/or `downstream reference` to satisfy you needs.

This simple example of double cleavage protocols is actually the large-deletion case. There are more cases in double cleavage protocols, say the upstream inversion, the downstream inversion, and the duplication. For more details, see https://doi.org/10.1016/j.molcel.2018.06.021.

To run the example, install rearr from conda.
```shell
$ conda install bioconda::rearr
```

Then execute
```shell
$ makeTarget="example.alg" \
fastqFiles="example.fq" \
markerIndices="example.fa" \
minScores=20 \
refFile="example.ref" \
runWorkFlow.sh
```

`runWorkFlow.sh` use `make` under the hood. Thus, the extension `.alg` is significant in `makeTarget="example.alg"`. `minScores=20` specify the necessary score for reads to hit a marker in `markerIndices="example.fa"`. You should adjust `minScores` according to the size of the markers. The score is calculated by the bowtie2 flags `--ma 1 --mp 2,2 --rdg 3,1 --rfg 3,1`. For more details, see [`demultiplex.sh`][demultiplex.sh.md].

The output (i.e. `makeTarget`) `example.alg` is a file with three rows as a group. Each three-row group corresponds to the alignment results of a read in `example.fq`. Note that reads in `example.fq` not hitting any markers in `example.fa` are not mapped, thereby not appearing in `example.alg`. Each three-row group has the follow form.

```
┌──────────────────┬──────────────────┬───────────────────────────────────┬───────────────────────────────────────┬──────────────────────────────────────────────────────────────────────┬──────────────────────────────────────────────────────────────┬────────────────────────────────────────────────────────────────────┬────────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────────────────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────────┬────────────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────────────────────────────────────────┬──────────────────────────────────────────────────────────────┬─────────────────────────────────────────┬─────────────────────────────────────────┬──────────────────────────────────────────────────────────────────────────┐
| query read index | query read count | alignment score | reference index | upstream unaligned part of query read | start position of the upstream alignment block in upstream reference | start position of the upstream alignment block in query read | end position of the upstream alignment block in upstream reference | end position of the upstream alignment block in query read | unaligned part of query read between the upstream and downstream alignment blocks | start position of the downstream alignment block in downstream reference + length of upstream reference | start position of the downstream alignment block in query read | end position of the downstream alignment block in downstream reference + length of upstream reference | end position of the downstream alignment block in query read | downstream unaligned part of query read | cleavage position in upstream reference | cleavage position in downstream reference + length of upstream reference |
├──────────────────┴──────────────────┴───────────────────────────────────┴───────────────────────────────────────┴──────────────────────────────────────────────────────────────────────┴──────────────────────────────────────────────────────────────┴────────────────────────────────────────────────────────────────────┴────────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────────────────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────┴────────────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────────────────────────────────────────┴──────────────────────────────────────────────────────────────┴─────────────────────────────────────────┴─────────────────────────────────────────┴──────────────────────────────────────────────────────────────────────────┤
| reference alignment line                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
├──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
| query alignment line                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

The second row `reference alignment line` and the third row `query alignment line` stacked together give the full detail of the chimeric aligment from query read to the upstream and downstream references. The first row contains information of the chimeric alignment.

- Duplicated query reads in `example.fq` are removed. The resulting unique query reads are listed by dictionary order in `example.noDup` (the name is determined from `example.alg` with extension `.alg` replaced by `.noDup`) together with their duplicating counts. `query read index` is the 1-based row number of the unique query read in `example.noDup`. `query read count` is the duplicating count in `example.fq` of the unique query read.
- `alignment score` is the score of the chimeric alignment (the higher the better).
- `reference index` is the 0-based row number of the reference in `example.ref`.
- `upstream unaligned part of query read` is the part at the start of the query which cannot be mapped to the upstream reference. This usually happens when the upstream reference does not cover the upstream end of the query read (say truncated) or the query read contains unalignable parts (say adapter).
- `start position of the upstream alignment block in upstream reference` is understood literally.
- `start position of the upstream alignment block in query read` is understood literally.
- `end position of the upstream alignment block in upstream reference` is understood literally.
- `end position of the upstream alignment block in query read` is understood literally.
- `unaligned part of query read between the upstream and downstream alignment blocks` is the random (or unpredictable) insertion in acadamic articles.
- `start position of the downstream alignment block in downstream reference + length of upstream reference` is start position of the downstream alignment block in the concatenate of the upstream and downstream references.
- `start position of the downstream alignment block in query read` is understood literally.
- `end position of the downstream alignment block in downstream reference + length of upstream reference` is end position of the downstream alignment block in the concatenate of the upstream and downstream references.
- `end position of the downstream alignment block in query read` is understood literally.
- `downstream unaligned part of query read` similar as `upstream unaligned part of query read` but for the unaligned part at the end of the query reads.
- `cleavage position in upstream reference` is understood literally.
- `cleavage position in downstream reference + length of upstream reference` is cleavage position in the concatenate of the upstream and downstream references.

[rearr]: https://github.com/ljw20180420/rearr
[example.fa]: ./assets/example/example.fa
[example.ref]: ./assets/example/example.ref
[example.fq]: ./assets/example/example.fq
[demultiplex.sh.md]: ./auxilary/demultiplex.sh.html
