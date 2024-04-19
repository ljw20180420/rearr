#!/bin/bash
# Usage: getSxCsvFileRef.sh csvfile genomeref [ext1up ext1down ext2up ext2down]
csvfile=$1
getfastagenome=$2
bowtie2genome=${getfastagenome%.*}
ext1up=${3:-50}
ext1down=${4:-0}
ext2up=${5:-10}
ext2down=${6:-100}

# For NAA csvfiles, 17~18bp of target is "TT", which should be replaced by "CC" in order to map genome. After find the genome location, the retrieved reference need replace "GG" (target and reference always have opposite strands, so "CC" becomes "GG") back to "AA".
getSxCsvFileTarget.pl "${csvfile}" | bowtie2 --quiet --mm -x "${bowtie2genome}" -r -U - 2> /dev/null | samtools view | gawk -f sxTargetSam2Bed.awk -- ${ext1up} ${ext1down} ${ext2up} ${ext2down} | bedtools getfasta -s -fi "${getfastagenome}" -bed - | sed '1~2d' | getSxRefFile.pl ${ext1up} ${ext2up} ${csvfile: -6:1}
