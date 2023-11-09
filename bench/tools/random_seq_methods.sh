#!/bin/bash
generate_random_DNA()
{
    # usage: generate_random_DNA length
    tr -dc ACGT </dev/urandom | head -c $1; echo ''
}


reflen=200
sglen=20
cut=$(expr $reflen / 2)
sgstart=$(expr $cut - $sglen + 3)
ext=$(expr $reflen / 20)
ref=$(generate_random_DNA $reflen)
ref=$(printf $ref | head -c$(($sgstart+$sglen+1)))GG$(printf $ref | tail -c$(($reflen-$sgstart-$sglen-3)))
sed '1i >ref' <<<$ref >bench/ref.fa
sgRNA=$(printf $ref | cut -c$(($sgstart+1))-$(($sgstart+$sglen)))

readnum=1000
bench/tools/random_seq_methods.py $ref 0.02 $readnum >bench/random.seq
cut -f1 bench/random.seq | perl -nE 'say "\@seq" . $. . "\n" . $_ . "+\n" . "~" x (length($_) - 1)' > bench/random.fq
cut -f2- bench/random.seq > bench/random.seq2; mv bench/random.seq2 bench/random.seq

rm -rf bench/rearr; mkdir -p bench/rearr
cp bench/random.fq bench/rearr/random.fq
rearr_run.sh bench/rearr/random.fq $ref $sgRNA $ext $ext
tail -n+2 bench/rearr/random.fq.table.$cut.$ext.$ext | cut -f5,7,15 | awk -v OFS="\t" -v cut=$cut -v ext=$ext '{$2 -= 2 * ext; print $1,$2,$3}' | paste <(cut -f1 bench/rearr/random.fq.count) - | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' bench/rearr/random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V > bench/rearr.arr

rm -rf bench/RESSO
CRISPResso -r1 bench/random.fq -a $ref -g $sgRNA -o bench/RESSO -amas 0
unzip $(find bench/RESSO -name "*.zip") -d bench/RESSO
bench/tools/ARRANGE_RESULTS.py <bench/RESSO/Alleles_frequency_table.txt CRISPResso | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' bench/random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V >bench/CRISPResso.arr

mkdir -p bench/CRVS
bwa index -p bench/CRVS/ref bench/ref.fa
bwa mem -v 3 -T -9999 bench/CRVS/ref bench/random.fq | samtools sort -o bench/CRVS/random.bam
samtools index -b bench/CRVS/random.bam
Rscript bench/tools/variants.r $ref bench/CRVS/random.bam $sgstart $sglen | sort -k1,1V | bench/tools/ARRANGE_RESULTS.py CrisprVariants $sgstart | awk -v OFS="\t" -v cut=$cut 'NF==1{print $0,cut,cut,0} NF>1{print}' >bench/CrisprVariants.arr

rm -rf bench/AMP; mkdir -p bench/AMP
primer=$(cut -c $(expr $reflen / 4 - 5)-$(expr $reflen / 4) <<<$ref)
sed "2~4{s/^/$primer/}; 4~4{s/^/~~~~~~/}" bench/random.fq >bench/AMP/random.fq.primer
amplicon=$(cut -c $(expr $reflen / 4 - 5)-$cut <<<$ref | tr '[:upper:]' '[:lower:]')$(cut -c $(expr $cut + 1) <<<$ref | tr '[:lower:]' '[:upper:]')$(cut -c $(expr $cut + 2)-$reflen <<<$ref | tr '[:upper:]' '[:lower:]')
printf "ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor\n" >bench/AMP/config.csv
printf "ID_1,barcode_1,random.fq.primer,,group_1,0,%s,%s,,0,%s," $sgRNA $primer $amplicon >>bench/AMP/config.csv
Rscript bench/tools/AMP.r bench/AMP/config.csv bench/AMP/ bench/AMP/
bench/tools/ARRANGE_RESULTS.py amplican $(expr $reflen / 4 - 6) <bench/AMP/alignments/alignments.txt | sort -k1,1 | join -t $'\t' -1 1 -2 1 <(sed -nr '1~4N; s/\n/\t/p' bench/random.fq | sort -k1,1) - | cut -f2- | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' bench/random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V >bench/amplican.arr

rm -rf bench/ADIV; mkdir -p bench/ADIV
sed '2~4{s/^/AAAAAATGTAAAACGACGGCCAGT/; s/$/aagacac/}; 4~4{s/^/~~~~~~~~~~~~~~~~~~~~~~~~/; s/$/~~~~~~~/}' bench/random.fq >bench/ADIV/random.for.fq
perl -pE 'if ($. % 4 == 2){tr/acgtACGT/TGCATGCA/} if ($. % 2 == 0){chomp; $_ = scalar reverse; $_ .= "\n"}' bench/ADIV/random.for.fq >bench/ADIV/random.rev.fq
novoindex bench/ADIV/ref.nix bench/ref.fa
novoalign -t 0,6 -o SAM -d bench/ADIV/ref.nix -f bench/ADIV/random.for.fq bench/ADIV/random.rev.fq | samtools view -o bench/ADIV/random.bam
printf "plate0\tA0\trefin_ref_45-200\tM\t0\tAAAAAA\n" >bench/ADIV/random.bar
rm -rf bench/ADIV/Workdir_target_ bench/ADIV/Output_target_
cd bench/ADIV; ../ampliconDIVider/ampliconDIVider_driver.sh -b random.bar -f ../ref.fa -x ref.nix random.bam; cd ../..
samtools view bench/ADIV/Workdir_target_/random.aligned.bam | bench/tools/ARRANGE_RESULTS.py "ampliconDIVider" | sort -k1,1V >bench/ampliconDIVider.arr

# rm -rf bench/ZhangFeng; mkdir -p bench/ZhangFeng
# printf "random\trandom.fq\t%s\t%s\n" $sgRNA $ref >bench/ZhangFeng/samplefile 
# ~/miniconda3/envs/py27/bin/python Screening_Protocols_manuscript/calculate_indel.py -i ZhangFeng/samplefile -no-m


for soft in $(printf "rearr CRISPResso CrisprVariants amplican ampliconDIVider")
do
    join -t $'\t' -1 1 -2 1 <(sort -k1,1 bench/random.seq) <(sort -k1,1 bench/$soft.arr) | sort -k1,1V | bench/tools/compare_indel.py $ref >bench/$soft.diff
done

bench/tools/draw_bench.py $readnum "bench/rearr.diff" "bench/CRISPResso.diff" "bench/CrisprVariants.diff" "bench/amplican.diff" "bench/ampliconDIVider.diff"