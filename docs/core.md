---
title: About
---

# `rearrangement`

`rearrangement` is the core chimeric alignment engine of [`rearr`][rearr].

{% highlight shell %}
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
{% endhighlight %}

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
                <td>ref</td>
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
- The second and third lines together form the actual alignment.
```

[rearr]: https://github.com/ljw20180420/rearr
