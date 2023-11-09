#!/bin/bash
generate_random_DNA()
{
    # usage: generate_random_DNA length
    tr -dc ACGT </dev/urandom | head -c $1; echo ''
}


reflen=200
sgstart=83
sglen=20
ref=$(generate_random_DNA $reflen)
ref=$(printf $ref | head -c$(($sgstart+$sglen+1)))GG$(printf $ref | tail -c$(($reflen-$sgstart-$sglen-3)))
sed '1i >ref' <<<$ref >ref.fa
sgRNA=$(printf $ref | cut -c$(($sgstart+1))-$(($sgstart+$sglen)))

readnum=1000
random_seq_methods.py $ref 0.02 $readnum > random.seq
cut -f1 random.seq | perl -nE 'say "\@seq" . $. . "\n" . $_ . "+\n" . "~" x (length($_) - 1)' > random.fq
cut -f2- random.seq > random.seq2; mv random.seq2 random.seq

rearr_run.sh random.fq $ref $sgRNA $(expr ${#ref} / 20) $(expr ${#ref} / 20)
join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' random.fq | sort -k2,2) <(cut -f1 random.fq.count | paste - <(tail -n+2 random.fq.table.$(expr ${#ref} / 2).$(expr ${#ref} / 20).$(expr ${#ref} / 20) | cut -f5,7,15 | awk -v OFS="\t" -v cut=$(expr ${#ref} / 2) -v ext1=$(expr ${#ref} / 20) -v ext2=$(expr ${#ref} / 20) '{$2 -= ext1 + ext2; print $1,$2,$3}') | sort -k1,1) | cut -f2- | sort -k1,1V > rearr.arr

rm -rf RESSO
CRISPResso -r1 random.fq -a $ref -g $sgRNA -o RESSO -amas 0
unzip $(find RESSO -name "*.zip") -d RESSO
join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' random.fq | sort -k2,2) <(ARRANGE_RESULTS.py <RESSO/Alleles_frequency_table.txt CRISPResso | sort -k1,1) | cut -f2- | sort -k1,1V >CRISPResso.arr

mkdir -p CRVS
bwa index -p CRVS/ref ref.fa
bwa mem -v 3 -T -9999 CRVS/ref random.fq | samtools sort -o CRVS/random.bam
samtools index -b CRVS/random.bam
Rscript variants.r $ref CRVS/random.bam $sgstart $sglen | sort -k1,1V > CRVS/random.crvs
ARRANGE_RESULTS.py CrisprVariants $sgstart <CRVS/random.crvs | awk -v OFS="\t" -v cut=$(expr ${#ref} / 2) 'NF==1{print $0,cut,cut,0} NF>1{print}' >CrisprVariants.arr

rm -rf AMP; mkdir -p AMP
primer=$(cut -c $(expr $reflen / 4 - 5)-$(expr $reflen / 4) <<<$ref)
sed "2~4{s/^/$primer/}; 4~4{s/^/~~~~~~/}" random.fq >AMP/random.fq.primer
amplicon=$(cut -c 45-$(expr $reflen / 2) <<<$ref | tr '[:upper:]' '[:lower:]')$(cut -c $(expr $reflen / 2 + 1) <<<$ref | tr '[:lower:]' '[:upper:]')$(cut -c $(expr $reflen / 2 + 2)-$reflen <<<$ref | tr '[:upper:]' '[:lower:]')
printf "ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor\n" > AMP/config.csv
printf "ID_1,barcode_1,random.fq.primer,,group_1,0,%s,%s,,0,%s," $sgRNA $primer $amplicon >> AMP/config.csv
Rscript AMP.r AMP/config.csv AMP/ AMP/
ARRANGE_RESULTS.py amplican $(expr $reflen / 4 - 6) <AMP/alignments/alignments.txt | sort -k1,1 | join -t $'\t' -1 1 -2 1 <(sed -nr '1~4N; s/\n/\t/p' random.fq | sort -k1,1) - | cut -f2- | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V >amplican.arr

rm -rf ADIV; mkdir -p ADIV
sed '2~4{s/^/AAAAAATGTAAAACGACGGCCAGT/; s/$/aagacac/}; 4~4{s/^/~~~~~~~~~~~~~~~~~~~~~~~~/; s/$/~~~~~~~/}' random.fq >ADIV/random.for.fq
perl -pE 'if ($. % 4 == 2){tr/acgtACGT/TGCATGCA/} if ($. % 2 == 0){chomp; $_ = scalar reverse; $_ .= "\n"}' ADIV/random.for.fq >ADIV/random.rev.fq
novoindex ADIV/ref.nix ref.fa
novoalign -t 0,6 -o SAM -d ADIV/ref.nix -f ADIV/random.for.fq ADIV/random.rev.fq | samtools view -o ADIV/random.bam
printf "plate0\tA0\trefin_ref_45-200\tM\t0\tAAAAAA\n" > ADIV/random.bar
rm -rf Workdir_target_ Output_target_
ampliconDIVider/ampliconDIVider_driver.sh -b ADIV/random.bar -f ref.fa -x ADIV/ref.nix ADIV/random.bam
samtools view Workdir_target_/random.aligned.bam | ARRANGE_RESULTS.py "ampliconDIVider" | sort -k1,1V >ampliconDIVider.arr

# rm -rf ZhangFeng; mkdir -p ZhangFeng
# printf "random\trandom.fq\t%s\t%s\n" $sgRNA $ref > ZhangFeng/samplefile 
# ~/miniconda3/envs/py27/bin/python Screening_Protocols_manuscript/calculate_indel.py -i ZhangFeng/samplefile -no-m


for soft in $(printf "rearr CRISPResso CrisprVariants amplican ampliconDIVider")
do
    join -t $'\t' -1 1 -2 1 <(sort -k1,1 random.seq) <(sort -k1,1 $soft.arr) | sort -k1,1V | compare_indel.py $ref >$soft.diff
done

draw_bench.py $readnum "rearr.diff" "CRISPResso.diff" "CrisprVariants.diff" "amplican.diff" "ampliconDIVider.diff"