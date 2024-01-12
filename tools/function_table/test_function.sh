#!/bin/bash

path="tools/function_table/single.100.0.01"
ref="$(sed '1d' tools/function_table/single.100.0.01/ref.fa)"
cut=$((${#ref} / 2))
sglen=20
sgstart=$(expr $cut - $sglen + 3)
sgRNA=$(printf $ref | cut -c$(($sgstart+1))-$(($sgstart+$sglen)))

echo "RESSO"
rm -rf $path/RESSO
/home/ljw/py310/bin/CRISPResso -r1 "$path/random.fq.gz" -a "$ref" -g "$sgRNA" -o $path/RESSO -amas 0

echo "AMP"
rm -rf $path/AMP; mkdir -p $path/AMP
primer=$(cut -c $(expr ${#ref} / 4 - 5)-$(expr ${#ref} / 4) <<<$ref)
zcat $path/random.fq.gz | sed "2~4{s/^/$primer/}; 4~4{s/^/~~~~~~/}" >$path/AMP/random.fq.primer
gzip $path/AMP/random.fq.primer
amplicon=$(cut -c $(expr ${#ref} / 4 - 5)-$cut <<<$ref | tr '[:upper:]' '[:lower:]')$(cut -c $(expr $cut + 1) <<<$ref | tr '[:lower:]' '[:upper:]')$(cut -c $(expr $cut + 2)-${#ref} <<<$ref | tr '[:upper:]' '[:lower:]')
printf "ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor\n" >$path/AMP/config.csv
printf "ID_1,barcode_1,random.fq.primer.gz,,group_1,0,%s,%s,,0,%s," $sgRNA $primer $amplicon >>$path/AMP/config.csv
Rscript bench/tools/AMP.r $path/AMP/config.csv $path/AMP/ $path/AMP/

echo "CRVS"
mkdir -p $path/CRVS
bwa index -p $path/CRVS/ref $path/ref.fa
bwa mem -v 3 -T -9999 $path/CRVS/ref $path/random.fq.gz | samtools sort -o $path/CRVS/random.bam
samtools index -b $path/CRVS/random.bam
Rscript bench/tools/variants.r $ref $path/CRVS/random.bam $(expr $cut - ${#sgRNA} + 3) ${#sgRNA} >$path/CRVS/result

echo "ADIV"
rm -rf $path/ADIV; mkdir -p $path/ADIV
barcode="AAAAAATGTAAAACGACGGCCAGT"
zcat $path/random.fq.gz | sed '2~4{s/^/'"$barcode"'/; s/$/aagacac/}; 4~4{s/^/~~~~~~~~~~~~~~~~~~~~~~~~/; s/$/~~~~~~~/}' >$path/ADIV/random.for.fq
perl -pE 'if ($. % 4 == 2){tr/acgtACGT/TGCATGCA/} if ($. % 2 == 0){chomp; $_ = scalar reverse; $_ .= "\n"}' $path/ADIV/random.for.fq >$path/ADIV/random.rev.fq
# gzip $path/ADIV/random.for.fq $path/ADIV/random.rev.fq
printf "plate0\tA0\trefin_ref_45-200\tM\t0\tAAAAAA\n" >$path/ADIV/random.bar
cd $path/ADIV
ADIVexu=$(realpath --relative-to=. $OLDPWD/bench/ampliconDIVider)/ampliconDIVider_driver.sh
cd - >/dev/null
novoindex $path/ADIV/ref.nix $path/ref.fa
# novoalign -t 0,6 -o SAM -d $path/ADIV/ref.nix -f $path/ADIV/random.for.fq.gz $path/ADIV/random.rev.fq.gz | samtools view -o $path/ADIV/random.bam
novoalign -t 0,6 -o SAM -d $path/ADIV/ref.nix -f $path/ADIV/random.for.fq $path/ADIV/random.rev.fq | samtools view -o $path/ADIV/random.bam
rm -rf $path/ADIV/Workdir_target_ $path/ADIV/Output_target_
cd $path/ADIV
$ADIVexu -b random.bar -f ../ref.fa -x ref.nix random.bam
cd -

echo "CRGR"
rm -rf $path/CRGR; mkdir -p $path/CRGR
cp $path/random.fq.gz $path/CRGR/random.for.fq.gz
zcat $path/CRGR/random.for.fq.gz | perl -pE 'if ($. % 4 == 2){tr/acgtACGT/TGCATGCA/} if ($. % 2 == 0){chomp; $_ = scalar reverse; $_ .= "\n"}' >$path/CRGR/random.rev.fq
gzip $path/CRGR/random.rev.fq
cp $path/ref.fa $path/CRGR/ref.fa
bench/CRISPR-GRANT/bin/indel_analysis -1:$path/CRGR/random.for.fq.gz -2:$path/CRGR/random.rev.fq.gz -r:$path/CRGR/ref.fa -o:$path/CRGR/output

echo "ZhangFeng"
rm -rf $path/ZhangFeng; mkdir -p $path/ZhangFeng
ZhangFengRefStart=$(expr ${#ref} / 4)
zcat $path/random.fq.gz >$path/random.fq
printf "random,$path/random.fq,%s,%s\n" $sgRNA $(cut -c$(expr $ZhangFengRefStart + 1)-$(expr ${#ref} / 4 \* 3) <<<$ref) >$path/ZhangFeng/samplefile
bench/Screening_Protocols_manuscript/my_calculate_indel.py -i $path/ZhangFeng/samplefile -no-m -o $path/ZhangFeng/calc_indel_out.csv >$path/ZhangFeng/result

echo "SelfTarget"
rm -rf $path/SelfTargetResults; mkdir -p $path/SeltTargetResults
printf ">ref %d FORWARD\n%s\n" $(expr $cut + 3) $ref >$path/SeltTargetResults/oligo.fa
# bench/SelfTarget/indel_analysis/indelmap/indelmap $path/random.fq.gz $path/SeltTargetResults/oligo.fa $path/SeltTargetResults/outputfile 0
bench/SelfTarget/indel_analysis/indelmap/indelmap $path/random.fq $path/SeltTargetResults/oligo.fa $path/SeltTargetResults/outputfile 0