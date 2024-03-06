#!/bin/bash

# usage
# rearr_barcode_index.sh "genome_fasta"

get_reference()
{
    local csvfile=$1
    local getfastagenome=$2
    local bowtie2genome=${getfastagenome%.*}
    local ext1up=$3
    local ext1down=$4
    local ext2up=$5
    local ext2down=$6
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
}

getfastagenome=$1
ext1up=50
ext1down=10
ext2up=10
ext2down=100

csvpath="$(dirname $(realpath $0))/../csvfiles"
for csvfile in final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A2.csv final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A3.csv final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G2.csv final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G3.csv
do
    perl -anF, -E '$pribar=substr($F[1], length($F[1]) - 41, 39); $pribar=~tr/ACGT/TGCA/; say ">BC_" . $F[0], "\n", scalar reverse $pribar' "$csvpath/$csvfile" >"$csvpath/$csvfile.primer+barcode.fa"
    bowtie2-build -q "$csvpath/$csvfile.primer+barcode.fa" "$csvpath/$csvfile.primer+barcode"

    perl -anF, -E '$rev=scalar reverse $F[1]; $rev=~m/[acgt]/g; say ">BC_" . $F[0], "\n", substr($F[1], 0, length($F[1]) - pos($rev) + 1)' "$csvpath/$csvfile" >"$csvpath/$csvfile.sgRNA+scaffold.fa"
    bowtie2-build -q "$csvpath/$csvfile.sgRNA+scaffold.fa" "$csvpath/$csvfile.sgRNA+scaffold"

    get_reference "$csvpath/$csvfile" "$getfastagenome" $ext1up $ext1down $ext2up $ext2down >"$csvpath/$csvfile.ref12"
done