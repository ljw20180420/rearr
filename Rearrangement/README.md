# Basic Usage
```{bash}
rearrangement <inputFile 3<refFile [options]
suggest options: -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9
```

# Parameters
- -h/-help/--help: Display help.
- -s0: Mismatching score. (default: -3)
- -s1: Matching score for non-extension reference part. (default: +1)
- -s2: Matching score for extension reference part. (default: +1)
- -u: Gap-extending penalty. (default: -2)
- -v: Gap-opening penalty. (default: -5)
- -ru: Gap-extending penalty for unaligned reference end. (default: 0)
- -rv: Gap-opening penalty for unaligned reference end. (default: 0)
- -qu: Gap-extending penalty for unaligned query part. (default: 0)
- -qv: Gap-opening penalty for unaligned query part. (default: 0)

# inputFile format
```
query1<tab>count1<tab>ref_id1<newline>
query2<tab>count2<tab>ref_id2<newline>
query3<tab>count3<tab>ref_id3<newline>
etc
```

# refFile format
```json
upBound11<tab>ref1seg1<tab>downBound11<tab>upBound12<tab>ref1seg2<tab>downBound12<tab>upBound13<tab>ref1seg3<tab>downBound13<tab>etc<newline>
upBound21<tab>ref2seg1<tab>downBound21<tab>upBound22<tab>ref2seg2<tab>downBound22<tab>upBound23<tab>ref2seg3<tab>downBound23<tab>etc<newline>
upBound31<tab>ref3seg1<tab>downBound31<tab>upBound32<tab>ref3seg2<tab>downBound32<tab>upBound33<tab>ref3seg3<tab>downBound33<tab>etc<newline>
etc
```

# outputFile format
```
index<tab>count<tab>score<tab>refId<tab>5'QueryUnmap<tab>ref1MapStart<tab>query1MapStart<tab>ref1MapEnd<tab>query1MapEnd<tab>randomInsertion<tab>ref2MapStart<tab>query2MapStart<tab>ref2MapEnd<tab>query2MapEnd<tab>3'QueryUnmap<tab>cut1<tab>cut2
```