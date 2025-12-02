#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
sxExtractMarker.sh csvfile >marker1 3>marker2
```

# Introduction
This is an in-house script to extract `marker`s from the `csvfile` of sx and lcy. The composition of the input csvfile is 
```
adapter(20bp) + sgRNA(20bp) + scaffold(83/93bp) + target(44bp) + 3bp + RCbarcode(18bp) + RCprimer(21bp)
```
Two `marker`s are output. `marker1` is output from `fd 1` (`stdout`). `marker1` consists of `primer` and `barcode`. `R2` is aligned to `marker1`. Then 44bp `target` is 3bp downstream to the end of `barcode` in `R2`. `marker2` is output from `fd 3`. `marker2` consists of `adapter`, `sgRNA` and `scaffold`. `R1` is aligned to `marker2`.

# Source
~~~bash
getSxCsvFilePrimer()
{
    # Usage: getSxCsvFilePrimer <csvfile
    rev | sed -r 's/^\s+//' | cut -c1-21 | dd conv=ucase 2>/dev/null | tr 'ACGT' 'TGCA'
}

getSxCsvFileBarcode()
{
    # Usage: getSxCsvFileBarcode <csvfile
    rev | sed -r 's/^\s+//' | cut -c22-39 | dd conv=ucase 2>/dev/null | tr 'ACGT' 'TGCA'
}

getSxCsvFileAdapter()
{
    # Usage: getSxCsvFileAdapter <csvfile
    cut -d, -f2 | cut -c1-20 | dd conv=ucase 2>/dev/null
}

getSxCsvFilesgRNA()
{
    # Usage: getSxCsvFilesgRNA <csvfile
    cut -d, -f2 | cut -c21-40 | dd conv=ucase 2>/dev/null
}

getSxCsvFileScaffold()
{
    # Usage: getSxCsvFileScaffold <csvfile
    cut -d, -f2 | sed -r 's/^[ACGTN]+//; s/[ACGTN]+\s+$//' | dd conv=ucase 2>/dev/null
}

getSxFaHead()
{
    # Usage: getSxFaHead <csvfile
    awk '{print ">" NR - 1}'
}

csvfile=$1
paste -d "" <(getSxCsvFilePrimer <"${csvfile}") <(getSxCsvFileBarcode <"${csvfile}") | paste -d "\n" <(getSxFaHead <"${csvfile}") - >&1

paste -d "" <(getSxCsvFileAdapter <"${csvfile}") <(getSxCsvFilesgRNA <"${csvfile}") <(getSxCsvFileScaffold <"${csvfile}") | paste -d "\n" <(getSxFaHead <"${csvfile}") - >&3
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
