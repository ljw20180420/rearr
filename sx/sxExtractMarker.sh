#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
```bash
$ sxExtractMarker.sh plasmid_file \
    >marker1 \
    3>marker2
```

<pre class="mermaid">
---
title: sxExtractMarker.sh
---
flowchart TD
    PF[(
        <h1>plasmid_file</h1>
        <table>
            <tr>
                <th>adapter#0040;20bp#0041; + sgRNA#0040;20bp#0041; + scaffold#0040;83/93bp#0041; + query#0040;44bp#0041; + 3bp + RCbarcode#0040;18bp#0041; + RCprimer#0040;21bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )] --> SEM[sxExtractMarker.sh]
    SEM --> MK1[(
        <h1>stdout</h1>
        <table>
            <tr>
                <th>primer#0040;21bp#0041; + barcode#0040;18bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )]
    SEM --> MK2[(
        <h1>fd3</h1>
        <table>
            <tr>
                <th>adapter#0040;20bp#0041; + sgRNA#0040;20bp#0041; + scaffold#0040;83/93bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )]
</pre>

- Extract `markers` from the in-house `plasmid_file`.
- `marker1` (primer(21bp) + barcode(18bp)) is output from `fd1` (`stdout`). `R2` is aligned to `marker1`. The 44bp `query` is 3bp downstream to the end of `barcode` in `R2`.
- `marker2` adapter(20bp) + sgRNA(20bp) + scaffold(83/93bp) is output from `fd3`. `R1` is aligned to `marker2`.

# Source
~~~bash
getSxPlasmidFilePrimer()
{
    # Usage: getSxPlasmidFilePrimer <plasmid_file
    rev | sed -r 's/^\s+//' | cut -c1-21 | dd conv=ucase 2>/dev/null | tr 'ACGT' 'TGCA'
}

getSxPlasmidFileBarcode()
{
    # Usage: getSxPlasmidFileBarcode <plasmid_file
    rev | sed -r 's/^\s+//' | cut -c22-39 | dd conv=ucase 2>/dev/null | tr 'ACGT' 'TGCA'
}

getSxPlasmidFileAdapter()
{
    # Usage: getSxPlasmidFileAdapter <plasmid_file
    cut -d, -f2 | cut -c1-20 | dd conv=ucase 2>/dev/null
}

getSxPlasmidFilesgRNA()
{
    # Usage: getSxPlasmidFilesgRNA <plasmid_file
    cut -d, -f2 | cut -c21-40 | dd conv=ucase 2>/dev/null
}

getSxPlasmidFileScaffold()
{
    # Usage: getSxPlasmidFileScaffold <plasmid_file
    cut -d, -f2 | sed -r 's/^[ACGTN]+//; s/[ACGTN]+\s+$//' | dd conv=ucase 2>/dev/null
}

getSxFaHead()
{
    # Usage: getSxFaHead <plasmid_file
    awk '{print ">" NR - 1}'
}

plasmid_file=$1
paste -d "" <(getSxPlasmidFilePrimer <"${plasmid_file}") <(getSxPlasmidFileBarcode <"${plasmid_file}") | paste -d "\n" <(getSxFaHead <"${plasmid_file}") - >&1

paste -d "" <(getSxPlasmidFileAdapter <"${plasmid_file}") <(getSxPlasmidFilesgRNA <"${plasmid_file}") <(getSxPlasmidFileScaffold <"${plasmid_file}") | paste -d "\n" <(getSxFaHead <"${plasmid_file}") - >&3
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
