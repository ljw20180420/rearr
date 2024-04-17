# Basic Usage
```{bash}
rearrangement <input_file 3<reference_file [options]
```

# Parameters
```{list}
-h/-help/--help: Display help.
-s0: Mismatching score. (default: -3)
-s1: Matching score for non-extension reference part. (default: +1)
-s2: Matching score for extension reference part. (default: +1)
-u: Gap-extending penalty. (default: -2)
-v: Gap-opening penalty. (default: -5)
-ru: Gap-extending penalty for unaligned reference end. (default: 0)
-rv: Gap-opening penalty for unaligned reference end. (default: 0)
-qu: Gap-extending penalty for unaligned query part. (default: 0)
-qv: Gap-opening penalty for unaligned query part. (default: 0)
```
