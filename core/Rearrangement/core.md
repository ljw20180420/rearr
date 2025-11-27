# Basic Usage
```{bash}
rearrangement (--help|-help|-h)
rearrangement <inputFile 3<refFile [options]
```

# `inputFile`
The one-line format is
```
query<tab>count<tab>refId<newline>
```
`query` is the query sequence. `count` is the duplicate number of the `query`. `refId` indicates the 0-based line number of reference in `refFline`

# `refFile`
The one-line format is
```json
s1<tab>ref1<tab>e1<tab>s2<tab>ref2<tab>e2<tab>...<tab>sN<tab>refN<tab>eN<newline>
```
`sM` and `eM` are the central range of `ref1`. The regions upstream to `sM` or downstream to `eM` are extension regions. The match bonus for extension regions are less than those for central region.

# `outputFile`
The three-line format is
```
index<tab>count<tab>score<tab>refId<newline>
refLine<newline>
queryLine<newline>
```
`index` is the 1-based alignment index as 1, 2, 3 and so on. `count` is copied from `inputFile`. `score` is the alignment score. `refId` is copied from `inputFile`. `refLine` and `queryLine` form the actual alignment.
