#!/bin/bash

# Usage: sxExtractSpliter.sh csvfiles
# csvfile = adapter(20bp) + sgRNA(20bp) + scaffold(83/93bp) + target(44bp) + 3bp + RCbarcode(18bp) + RCprimer(21bp)

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

for csvfile in "$@"
do
    paste -d "" <(getSxCsvFilePrimer <"${csvfile}") <(getSxCsvFileBarcode <"${csvfile}") | paste -d "\n" <(getSxFaHead <"${csvfile}") - >"${csvfile}.primer+barcode.fa"

    paste -d "" <(getSxCsvFileAdapter <"${csvfile}") <(getSxCsvFilesgRNA <"${csvfile}") <(getSxCsvFileScaffold <"${csvfile}") | paste -d "\n" <(getSxFaHead <"${csvfile}") - >"${csvfile}.adapter+sgRNA+scaffold.fa"
done