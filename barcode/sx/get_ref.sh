#!/bin/bash
# Usage: get_ref.sh csvfile fasta_genome ext1up ext1down ext2up ext2down
csvfile=$1
getfastagenome=$2
bowtie2genome=${getfastagenome%.*}
ext1up=${3:-50}
ext1down=${4:-10}
ext2up=${5:-10}
ext2down=${6:-100}

# For NAA csvfiles, 17~18bp of target is "TT", which should be replaced by "CC" in order to map genome. After find the genome location, the retrieved reference need replace "GG" (target and reference always have opposite strands, so "CC" becomes "GG") back to "AA".
perl -anF, -E '
    $rev = scalar reverse $F[1];
    $rev=~m/[acgt]/g;
    $target=substr($F[1], length($F[1]) - pos($rev) + 1, 44);
    substr($target, 16, 2)="CC";
    say $target;
' "$csvfile" | bowtie2 --quiet --mm -x "$bowtie2genome" -r -U - 2> /dev/null | samtools view | awk -F "\t" -v OFS="\t" -v ext1up="$ext1up" -v ext1down="$ext1down" -v ext2up="$ext2up" -v ext2down="$ext2down" '
    {
        qname = $1
        flag = $2
        if (int(flag / 16) % 2)
            refstrand = "+"
        else
            refstrand = "-"
        chr = $3
        pos = $4 - 1
        if (refstrand == "+")
        {
            target_length = substr($6, 1, length($6) - 1)
            cut = pos + target_length - 16 - 6
            print chr, cut - ext1up, cut + ext1down, qname, ".", "+"
            print chr, cut - ext2up, cut + ext2down, qname, ".", "+"
        }
        else
        {
            cut = pos + 16 + 6
            print chr, cut - ext1down, cut + ext1up, qname, ".", "-"
            print chr, cut - ext2down, cut + ext2up, qname, ".", "-"
        }
    }
' | bedtools getfasta -s -fi "$getfastagenome" -bed - | sed '1~2d' | perl -snE '
    if ($AG eq "A"){
        if ($. % 2){
            $spos = $ext1up + 4
        }
        else{
            $spos = $ext2up + 4
        }
        substr($_, $spos, 2) = "AA"
    }
    if ($. % 2){
        chomp;
        print $_
    }
    else{
        print "\t" . $_
    }
' -- -ext1up="$ext1up" -ext2up="$ext2up" -AG="${csvfile: -6:1}"