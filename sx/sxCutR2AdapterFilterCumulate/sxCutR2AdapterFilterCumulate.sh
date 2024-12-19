#!/bin/bash
# Remove spliterToMap from 5' and adapter from 3' of toMap. The remain must contain at least minToMapShear bases. After shearing, adjacent lines may be duplicate again. Use sxCumulateToMapCutAdaptSpliter.awk to collapse these duplicates.
# Usage: sxCutR2AdapterFilterCumulate.sh demultiplexFile minToMapShear >toMapFile

cutadaptPlain()
{
    # Usage: cutadaptPlain <plainseq 3'adapter
    # cutadapt does not accept plainseq. This function transform plainseq to fasta before feed to cutadapt, and then transform the fasta output back to plainseq
    # Input: plainseq
    # Output: 3' trimmed plainseq
    sed '=' | sed '1~2s/^/>s/' | cutadapt -a $1 - 2> /dev/null | sed '1~2d'
}

rmDupFile=$1
minToMapShear=$2
cut -f1 $rmDupFile | cutadaptPlain GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC | paste - <(cut -f3-4,8 $rmDupFile) | gawk -F "\t" -v OFS="\t" -v minToMapShear=$minToMapShear '
{
    if ($4 + minToMapShear <= length($1)) {
        print substr($1, $4 + 1), $2, $3
    }
}' | gawk -f sxCumulateToMapCutAdaptSpliter.awk
