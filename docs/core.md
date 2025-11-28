---
title: About
---

# `rearrangement`

`rearrangement` is the core chimeric alignment engine of [`rearr`][rearr].

```shell
$ rearrangement -h
### Basic Usage
rearrangement <input_file 3<reference_file

### Parameters
-h, -help, --help: Display help.
# Aligning Parameters
-s0: Mismatching score. (default: -6)
-s1: Matching score for non-extension reference part. (default: 4)
-s2: Matching score for extension reference part. (default: 2)
-u: Gap-extending penalty. (default: -3)
-v: Gap-opening penalty. (default: -9)
-ru: Gap-extending penalty for unaligned reference ends. (default: 0)
-rv: Gap-opening penalty for unaligned reference ends. (default: 0)
-qu: Gap-extending penalty for unaligned query parts. (default: 0)
-qv: Gap-opening penalty for unaligned query parts. (default: -5)
```

<pre class="mermaid">
---
title: rearrangement
---
flowchart TD
    QUERY[(
        <h1>input_file</h1>
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
    )] --> REARR[rearrangement]
    REF[(
        <h1>reference_file</h1>
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
    )] --> REARR
    REARR --> ALG[(
        <h1>stdout</h1>
        <table>
            <tr>
                <td>idx</td>
                <td>#</td>
                <td>score</td>
                <td>id</td>
            </tr>
            <tr>
                <td>ref1</td>
                <td>ref2</td>
            </tr>
            <tr>
                <td>query</td>
            </tr>
        </table>
    )]
</pre>

## `input_file`

- `query` is the query sequence.
- `#` is the duplication number of the `query`.
- `id` is the 0-based line number of reference in `reference_file`.

## `reference_file`

- [`rearr`][rearr] takes two references (`ref1` and `ref2`) as input. This is useful in analyzing the [double cleavage CRISPR experiments](https://doi.org/10.1016/j.molcel.2018.06.021), say to delete a large portion of the genome.
- The upstream part of query is aligned to`ref1`, which in the case of two-cleavage deletion, is the sequence around the upstream cleavage site, The actual `ref1` may depend on the repair junction (deletion, inversion or duplication).
- The downstream part of query is aligned to`ref2`, which in the case of two-cleavage deletion, is the sequence around the downstream cleavage site.
- For the single cleavage case, `ref2` just repeats `ref1`.
- The region between `start1` and `end1` is the non-extension region of `ref1`. The regions upstream to `start1` or downstream to `end1` are extension regions of `ref1`. The matching score (`s2`) for extension regions is lower than that (`s1`) for the non-extension region.
- The region between `start2` and `end2` is the non-extension region of `ref2`. The regions upstream to `start2` or downstream to `end2` are extension regions of `ref2`.
- The reason to distinguish extension and non-extension regions is to avoid dummy templated insertion induced by small mutations away from the cleavage site for the single cleavage case, as shown in the following video.

<video width="854" height="480" controls>
  <source src="/assets/videos/MainScene_0007_ScoreCrossCleavage.mp4" type="video/mp4">
</video>

## `stdout`

Every three lines of the standard output represents a single alignment.

- The first line is a header line.
  - `idx` is the 1-based line number of `query` in `input_file`.
  - `#` is the duplication number of the `query`.
  - `score` is the alignment score. 
  - `id` is the 0-based line number of reference in `reference_file`.
- The second line is the sequence of the reference.
- The third line is the query with `idx`.
- The second and third lines together form the actual alignment, as shown in the following example.

```plaintext
1       1       157     9300
---aGTTGGCTAGTCAATACCTGAAGAGAGATTGGCCTGGAGTAAAAGC-TGAtaAAAGCTGATGATCGGAATGATTACAGGTAAATTAGTAGTTTTTGCCTATTTTCTTTAGAAACGGTTTTACTTAAAGCTATGTTACATATAGATAATGTAACACTCTAGt-------
CTG----------------------------TTGGCCTGGAGTAAAAGCATGAT----------GATCGGAATGATTACAGGTAAA------------------------------------------------------------------------------CAAAAAA
```

# `correct_micro_homology.awk`

Microhomology is common in CRISPR editing output. When microhomology happens, `rearrangement` cannot determine how to align `query` to `ref1` and `ref2`, as show in the following video.

<video width="854" height="480" controls>
  <source src="/assets/videos/MainScene_0008_MicroHomology.mp4" type="video/mp4">
</video>

`correct_micro_homology.awk` allows one to specify which end of the double strand break should be corrected toward the cleavage site up to the microhomology equivalence.

```shell
$ gawk -f correct_micro_homology.awk -- \
    correct_micro_homology/reference_file \
    correct_micro_homology/direction_file \
    < correct_micro_homology/rearrangement_file
```

<pre class="mermaid">
---
title: correct_micro_homology.awk
---
flowchart TD
    REF[(
        <h1>reference_file</h1>
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
    )] --> CMH[correct_micro_homology.awk]
    DIRECTION[(
        <h1>direction_file</h1>
        <table>
            <tr>
                <th>up/down</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )] --> CMH
    ALG[(
        <h1>rearrangement_file</h1>
        <table>
            <tr>
                <td>idx</td>
                <td>#</td>
                <td>score</td>
                <td>id</td>
            </tr>
            <tr>
                <td>ref1</td>
                <td>ref2</td>
            </tr>
            <tr>
                <td>query</td>
            </tr>
        </table>
    )] --> CMH
    CMH --> CORRECTED[(
        <h1>stdout</h1>
        <table>
            <tr>
                <td>idx</td>
                <td>#</td>
                <td>score</td>
                <td>id</td>
                <td>udangle</td>
                <td>rstart1</td>
                <td>qstart1</td>
                <td>rend1</td>
                <td>qend1</td>
                <td>random</td>
                <td>rstart2</td>
                <td>qstart2</td>
                <td>rend2</td>
                <td>qend2</td>
                <td>ddangle</td>
                <td>cut1</td>
                <td>ref1+cut2</td>
            </tr>
            <tr>
                <td>ref1</td>
                <td>ref2</td>
            </tr>
            <tr>
                <td>query</td>
            </tr>
        </table>
    )]
</pre>

- `reference_file` is the same as that takes by `rearrangement`.
- Each line in `direction_file` corresponds to each line in `reference_file`, which contains an `up` or `down` string to specify whether the upstream DSB end or the downstream DSB end should be corrected towards the cleavage site.
- `rearrangement_file` is the `stdout` of `rearrangement`.
- `stdout` of `correct_micro_homology.awk` is similar to `stdout` of `rearrangement` but with an extended header line.
  - `udangle` is the upstream unaligned part of `query`.
  - `rstart1` and `rend1` specifies the `ref1` range for the upstream block of the chimeric alignment.
  - `qstart1` and `qend1` specifies the `query` range for the upstream block of the chimeric alignment.
  - `random` is the unaligned part of `query` between the upstream and downstream block of the chimeric alignment..
  - `rstart2` and `rend2` specifies the `ref2` range for the downstream block pf the chimeric alignment.
  - `qstart2` and `qend2` specifies the `query` range for the downstream block of the chimeric alignment.
  - `ddangle` is the downstream unaligned part of `query`.

# Core

`rearrangement` and `correct_micro_homology.awk` forms the core part of [`rearr`][rearr]. They are generally piped together.
```shell
$ rearrangement \
    < input_file \
    3< reference_file |
  gawk -f correct_micro_homology.awk -- \
    reference_file \
    direction_file
```

# More than two blocks

`rearrangement` and `correct_micro_homology.awk` supports chimeric alignments with more than two blocks. The core part of [`rearr`][rearr] for multiple blocks is as follows.

<pre class="mermaid">
---
title: core
---
flowchart TD
    QUERY[(
        <h1>input_file</h1>
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
    )] --> REARR[rearrangement]
    REF[(
        <h1>reference_file</h1>
        <table>
            <tr>
                <th>start1</th>
                <th>ref1</th>
                <th>end1</th>
                <th>start2</th>
                <th>ref2</th>
                <th>end2</th>
                <th>...</th>
                <th>startN</th>
                <th>refN</th>
                <th>endN</th>
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
            </tr>
        </table>
    )] --> REARR

    REF --> CMH[correct_micro_homology.awk]
    DIRECTION[(
        <h1>direction_file</h1>
        <table>
            <tr>
                <th>up/down:1:2</th>
                <th>up/down:2:3</th>
                <th>...</th>
                <th>up/down:N-1:N</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )] --> CMH
    REARR --> CMH
    CMH --> CORRECTED[(
        <h1>stdout</h1>
        <table>
            <tr>
                <td>idx</td>
                <td>#</td>
                <td>score</td>
                <td>id</td>
                <td>dangle0</td>
                <td>rstart1</td>
                <td>qstart1</td>
                <td>rend1</td>
                <td>qend1</td>
                <td>dangle1</td>
                <td>rstart2</td>
                <td>qstart2</td>
                <td>rend2</td>
                <td>qend2</td>
                <td>dangle2</td>
                <td>...</td>
                <td>rstartN</td>
                <td>qstartN</td>
                <td>rendN</td>
                <td>qendN</td>
                <td>dangleN</td>
            </tr>
            <tr>
                <td>ref1</td>
                <td>ref2</td>
            </tr>
            <tr>
                <td>query</td>
            </tr>
        </table>
    )]
</pre>

Each row of `direction_file` has multiple fields corresponding to the junctions of adjacent references. The extended header of `stdout` of `correct_micro_homology.awk` contains information for all alignment blocks and all unaligned parts of `query`.

[rearr]: https://github.com/ljw20180420/rearr
