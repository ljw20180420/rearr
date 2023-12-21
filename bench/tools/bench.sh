#!/bin/bash

mode=${1:-"single"}
reflen=${2:-100}
probability=${3:-0.01}
readnum=${4:-10000}

generate_random_DNA()
{
    # usage: generate_random_DNA length
    tr -dc ACGT </dev/urandom | head -c $1; echo ''
}

generate_reference()
{
    local reflen=$1
    local sglen=$2
    local cut=$(expr $reflen / 2)
    local sgstart=$(expr $cut - $sglen + 3)
    local ref=$(generate_random_DNA $reflen)
    local ref=$(printf $ref | head -c$(($sgstart+$sglen+1)))GG$(printf $ref | tail -c$(($reflen-$sgstart-$sglen-3)))
    local sgRNA=$(printf $ref | cut -c$(($sgstart+1))-$(($sgstart+$sglen)))
    echo $cut $ref $sgRNA
}

get_time_memory()
{
    read usertime systime realtime memory <<<$(eval /bin/time $1 2>&1 | tee -a $path/bench.log | sed -nr 's/^([.0-9]+)user ([.0-9]+)system ([.:0-9]+)elapsed.*\(.* ([0-9]+)maxresident\)k$/\1 \2 \3 \4/p')
    realtime=$(awk -F ":" '{print $1 * 60 + $2}' <<<$realtime)
}

format_output()
{
    sed "s/^/$mode\t$reflen\t$probability\t$readnum\t$1\t$usertime\t$systime\t$realtime\t$memory\t/" | bench/tools/compare_indel.py $ref1 $ref2 $path/random.seq
}

if [ $mode != "single"  ] && [ $mode != "double" ]
then
    echo "mode is neither single nor double"
    exit 1
fi

if [ $mode = "double" ]
then
    path="bench/double"
    read cut1 ref1 sgRNA1 <<<$(generate_reference $reflen 20)
    read cut2 ref2 sgRNA2 <<<$(generate_reference $reflen 20)
    ref=$(printf $ref1 | head -c$cut1)$(printf $ref2 | tail -c$((${#ref2} - $cut2)))
    cut=$(expr ${#ref} / 2)
    sgRNA=$(printf $ref | cut -c$(expr $cut - ${#sgRNA1} + 4)-$(expr $cut + 3))
else
    path="bench/single"
    read cut ref sgRNA <<<$(generate_reference $reflen 20)
    ref1=$ref
    cut1=$cut
    sgRNA1=$sgRNA
    ref2=$ref
    cut2=$cut
    sgRNA2=$sgRNA
fi
>$path/bench.log
sed '1i >ref' <<<$ref >$path/ref.fa
sed '1i >ref1' <<<$ref1 >$path/ref12.fa
sed '1i >ref2' <<<$ref2 >>$path/ref12.fa

bench/tools/random_seq_methods.py $ref1 $ref2 $probability $readnum >$path/random.seq
cut -f1 $path/random.seq | perl -nE 'say "\@seq" . $. . "\n" . $_ . "+\n" . "~" x (length($_) - 1)' >$path/random.fq
cut -f2- $path/random.seq > $path/random.seq2; mv $path/random.seq2 $path/random.seq
sed "s/^/$mode\t$reflen\t$probability\t$readnum\tsimulate\t*\t*\t*\t*\t/" $path/random.seq >&3

rm -rf $path/rearr; mkdir -p $path/rearr
cp $path/random.fq $path/rearr/random.fq
get_time_memory "rearr_run.sh $path/rearr/random.fq $ref1 $ref2 $sgRNA1 $sgRNA2"
tail -n+2 $path/rearr/random.fq.table.$cut1.$cut2.${#ref1} | awk -F "\t" -v OFS="\t" -v mid=$((${#ref1} - $cut1 + $cut2)) '{$10 -= mid; print $7,$10,$8, $11}' | paste <(cut -f1 $path/rearr/random.fq.count) - | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' $path/rearr/random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V | format_output rearr >&3

rm -rf $path/RESSO
get_time_memory "/home/ljw/py310/bin/CRISPResso -r1 $path/random.fq -a $ref -g $sgRNA -o $path/RESSO -amas 0"
unzip -q -o $(find $path/RESSO -name "*.zip") -d $path/RESSO
bench/tools/ARRANGE_RESULTS.py <$path/RESSO/Alleles_frequency_table.txt CRISPResso $cut | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' $path/random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V | format_output RESSO >&3

rm -rf $path/AMP; mkdir -p $path/AMP
primer=$(cut -c $(expr ${#ref} / 4 - 5)-$(expr ${#ref} / 4) <<<$ref)
sed "2~4{s/^/$primer/}; 4~4{s/^/~~~~~~/}" $path/random.fq >$path/AMP/random.fq.primer
amplicon=$(cut -c $(expr ${#ref} / 4 - 5)-$cut <<<$ref | tr '[:upper:]' '[:lower:]')$(cut -c $(expr $cut + 1) <<<$ref | tr '[:lower:]' '[:upper:]')$(cut -c $(expr $cut + 2)-${#ref} <<<$ref | tr '[:upper:]' '[:lower:]')
printf "ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor\n" >$path/AMP/config.csv
printf "ID_1,barcode_1,random.fq.primer,,group_1,0,%s,%s,,0,%s," $sgRNA $primer $amplicon >>$path/AMP/config.csv
get_time_memory "Rscript bench/tools/AMP.r $path/AMP/config.csv $path/AMP/ $path/AMP/"
bench/tools/ARRANGE_RESULTS.py amplican $cut $(expr ${#ref} / 4 - 6) ${#primer} <$path/AMP/alignments/alignments.txt | paste <(sed -n '2~4p' $path/AMP/random.fq.primer | sort | uniq -c | sort -s -k1,1nr | awk '{print $2}') - | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4s/^@//; N; s/\n/\t/p' $path/AMP/random.fq.primer | sort -k2,2) - | cut -f2- | sort -k1,1V | format_output AMP >&3

mkdir -p $path/CRVS
get_time_memory "bash -c \"bwa index -p $path/CRVS/ref $path/ref.fa; bwa mem -v 3 -T -9999 $path/CRVS/ref $path/random.fq | samtools sort -o $path/CRVS/random.bam; samtools index -b $path/CRVS/random.bam; Rscript bench/tools/variants.r $ref $path/CRVS/random.bam $(expr $cut - ${#sgRNA1} + 3) ${#sgRNA} >$path/CRVS/result\""
cut -f1 $path/CRVS/result | sort | join -t $'\t' -1 1 -2 1 - <(samtools view $path/CRVS/random.bam | sort -k1,1) | bench/tools/ARRANGE_RESULTS.py CrisprVariants $cut | sort -k1,1V | format_output CRVS >&3

rm -rf $path/ADIV; mkdir -p $path/ADIV
barcode="AAAAAATGTAAAACGACGGCCAGT"
sed '2~4{s/^/'"$barcode"'/; s/$/aagacac/}; 4~4{s/^/~~~~~~~~~~~~~~~~~~~~~~~~/; s/$/~~~~~~~/}' $path/random.fq >$path/ADIV/random.for.fq
perl -pE 'if ($. % 4 == 2){tr/acgtACGT/TGCATGCA/} if ($. % 2 == 0){chomp; $_ = scalar reverse; $_ .= "\n"}' $path/ADIV/random.for.fq >$path/ADIV/random.rev.fq
printf "plate0\tA0\trefin_ref_45-200\tM\t0\tAAAAAA\n" >$path/ADIV/random.bar
read usertime systime realtime memory <<<$(/bin/time bash -c 'novoindex '$path'/ADIV/ref.nix '$path'/ref.fa; novoalign -t 0,6 -o SAM -d '$path'/ADIV/ref.nix -f '$path'/ADIV/random.for.fq '$path'/ADIV/random.rev.fq | samtools view -o '$path'/ADIV/random.bam; rm -rf '$path'/ADIV/Workdir_target_ '$path'/ADIV/Output_target_; cd '$path'/ADIV; $(realpath --relative-to=. $OLDPWD/bench/ampliconDIVider)/ampliconDIVider_driver.sh -b random.bar -f ../ref.fa -x ref.nix random.bam; cd -' 2>&1 | sed -nr 's/^([.0-9]+)user ([.0-9]+)system ([.:0-9]+)elapsed.*\(.* ([0-9]+)maxresident\)k$/\1 \2 \3 \4/p')
realtime=$(awk -F ":" '{print $1 * 60 + $2}' <<<$realtime)
samtools view $path/ADIV/Workdir_target_/random.aligned.bam | bench/tools/ARRANGE_RESULTS.py "ampliconDIVider" $cut | sort -k1,1V | format_output ADIV >&3

rm -rf $path/CRGR; mkdir -p $path/CRGR
cp $path/random.fq $path/CRGR/random.for.fq
perl -pE 'if ($. % 4 == 2){tr/acgtACGT/TGCATGCA/} if ($. % 2 == 0){chomp; $_ = scalar reverse; $_ .= "\n"}' $path/CRGR/random.for.fq >$path/CRGR/random.rev.fq
cp $path/ref.fa $path/CRGR/ref.fa
read usertime systime realtime memory <<<$(/bin/time bench/CRISPR-GRANT/bin/indel_analysis -1:$path/CRGR/random.for.fq -2:$path/CRGR/random.rev.fq -r:$path/CRGR/ref.fa -o:$path/CRGR/output 2>&1 | sed -nr 's/^([.0-9]+)user ([.0-9]+)system ([.:0-9]+)elapsed.*\(.* ([0-9]+)maxresident\)k$/\1 \2 \3 \4/p')
realtime=$(awk -F ":" '{print $1 * 60 + $2}' <<<$realtime)
samtools view $path/CRGR/output/4.Mapping_sorted.bam | bench/tools/ARRANGE_RESULTS.py "CRISPR-GRANT" $cut | sort -k1,1V | format_output CRGR >&3

rm -rf $path/ZhangFeng; mkdir -p $path/ZhangFeng
ZhangFengRefStart=$(expr ${#ref} / 4)
printf "random,$path/random.fq,%s,%s\n" $sgRNA $(cut -c$(expr $ZhangFengRefStart + 1)-$(expr ${#ref} / 4 \* 3) <<<$ref) >$path/ZhangFeng/samplefile
get_time_memory "bench/Screening_Protocols_manuscript/my_calculate_indel.py -i $path/ZhangFeng/samplefile -no-m -o $path/ZhangFeng/calc_indel_out.csv >$path/ZhangFeng/result"
tail -n+2 $path/ZhangFeng/result | bench/tools/ARRANGE_RESULTS.py ZhangFeng $cut $ZhangFengRefStart | format_output ZhangFeng >&3

rm -rf $path/SelfTargetResults; mkdir -p $path/SeltTargetResults
printf ">ref %d FORWARD\n%s\n" $(expr $cut + 3) $ref >$path/SeltTargetResults/oligo.fa
get_time_memory "bench/SelfTarget/indel_analysis/indelmap/indelmap $path/random.fq $path/SeltTargetResults/oligo.fa $path/SeltTargetResults/outputfile 0"
tail -n+2 $path/SeltTargetResults/outputfile | bench/tools/ARRANGE_RESULTS.py SelfTarget $cut | format_output SelfTarget >&3

# cd $path/CRISPR_toolkit/Indel_searcher_2; ./Make_user_folder.sh
# cp ../../ref.fa Input/JaeWoo/Reference
# cp ../../random.fq Input/JaeWoo/FASTQ
# ./Run_cmd.sh
# cd ../../../

# rm -rf $path/canf; mkdir -p $path/canf
# sed -r -e '/^\s+outDir/c \        outDir = "output"' -e '/^\s+input/c \        input = "../random.fq"' -e '/^\s+r2file/c \        r2file = false' -e '/^\s+referenceFasta/c \        referenceFasta = "../ref.fa"' -e '/^\s+refOrganism/c \        refOrganism = [["ref", "other"]]' -e '/^\s+gRNAseq/c \        gRNAseq = [["ref", "'"$sgRNA"'"]]' -e '/^\s+selfRef/c \        selfRef = [["ref", false ]]' -e '/^\s+template_seq/c \        template_seq = "ref.fa"' -e '/^\s+protCutSite/c \        protCutSite = [["ref", -3]]' -e '/^\s+indexHuman/c \        indexHuman = "/home/ljw/hg19_with_bowtie2_index"' -e '/^\s+indexMouse/c \        indexMouse = "/home/ljw/mm9_with_bowtie2_index"' $path/crispr-a_nextflow/nf_input-example.config >$path/canf/random.config
# cd $path/canf; ../crispr-a_nextflow/nextflow-21.10.6-all run ../crispr-a_nextflow/crispr-a.nf -c random.config; cd ../..

# if [ $mode = "double" ]
# then
#     for soft in $(printf "rearr CRISPResso CrisprVariants amplican ZhangFeng SelfTarget")
#     do
#         join -t $'\t' -1 1 -2 1 <(sort -k1,1 $path/random.seq) <(sort -k1,1 $path/$soft.arr) | sort -k1,1V | bench/tools/compare_indel.py $ref >$path/$soft.diff
#     done

#     bench/tools/draw_bench.py $readnum "$path/rearr.diff" "$path/CRISPResso.diff" "$path/CrisprVariants.diff" "$path/amplican.diff" "$path/ZhangFeng.diff" "$path/SelfTarget.diff"
# else
#     for soft in $(printf "rearr CRISPResso CrisprVariants amplican ampliconDIVider CRISPR-GRANT ZhangFeng SelfTarget")
#     do
#         join -t $'\t' -1 1 -2 1 <(sort -k1,1 $path/random.seq) <(sort -k1,1 $path/$soft.arr) | sort -k1,1V | bench/tools/compare_indel.py $ref >$path/$soft.diff
#     done

#     bench/tools/draw_bench.py $readnum "$path/rearr.diff" "$path/CRISPResso.diff" "$path/CrisprVariants.diff" "$path/amplican.diff" "$path/ampliconDIVider.diff" "$path/CRISPR-GRANT.diff" "$path/ZhangFeng.diff" "$path/SelfTarget.diff"
# fi


