#!/bin/bash

get_reference()
{
    local csvfile=$1
    local bowtie2genome=$2
    local getfastagenome=$3
    local ext1up=$4
    local ext1down=$5
    local ext2up=$6
    local ext2down=$7
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

bowtie2genome=$1
getfastagenome=$2
ext1up=50
ext1down=10
ext2up=10
ext2down=100

for csvfile in final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A2.csv final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A3.csv final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G2.csv final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G3.csv
do
    perl -anF, -E '$pribar=substr($F[1], length($F[1]) - 41, 39); $pribar=~tr/ACGT/TGCA/; say ">BC_" . $F[0], "\n", scalar reverse $pribar' barcode/csvfiles/$csvfile >"barcode/csvfiles/$csvfile.primer+barcode.fa"
    bowtie2-build -q "barcode/csvfiles/$csvfile.primer+barcode.fa" "barcode/csvfiles/$csvfile.primer+barcode"

    perl -anF, -E '$rev=scalar reverse $F[1]; $rev=~m/[acgt]/g; say ">BC_" . $F[0], "\n", substr($F[1], 0, length($F[1]) - pos($rev) + 1)' barcode/csvfiles/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv >"barcode/csvfiles/$csvfile.sgRNA+scaffold.fa"
    bowtie2-build -q "barcode/csvfiles/$csvfile.sgRNA+scaffold.fa" "barcode/csvfiles/$csvfile.sgRNA+scaffold"

    get_reference "barcode/csvfiles/$csvfile" "$bowtie2genome" "$getfastagenome" $ext1up $ext1down $ext2up $ext2down >"barcode/csvfiles/$csvfile.ref12"
done