#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

title() {
    printf "\n----------\n%s\n----------\n" $1 >&2
}


to_fastq() {
    awk '
        {
            quality = sprintf("%*s", length($1), "")
            gsub(".", "~", quality)
            printf("@seq%d\n%s\n+\n%s\n", NR, $1, quality)
        }
    '
}

get_ref() {
    awk '
        {
            printf("%s%s\n", substr($2, 1, $3), substr($5, $4 + 1))
        }
    '
}

get_sgRNA() {
    awk '
        {
            printf("%s%s\n", substr($2, $3 - 16, 17), substr($5, $4 + 1, 3))
        }
    '
}

get_amplicon() {
    awk '
        {
            printf("%s%s%s\n", tolower(substr($2, $1 - 5, $3 - $1 + 6)), toupper(substr($2, $3 + 1, 1)), tolower(substr($5, $4 + 2)))
        }
    '
}

get_ref_ZhangFeng() {
    awk '
        {
            printf("%s%s\n", substr($2, $1 + 1, $3 - $1), substr($5, $4 + 1, $6 - $4))
        }
    '
}

row_print() {
    local software=$1
    local realtime=$2
    local usertime=$3
    local systime=$4
    local memory=$5
    sed "s/^/$mode\t$ref_len\t$probability\t$software\t$realtime\t$usertime\t$systime\t$memory\t/"
}

parse_rqpos() {
    awk '
        {
            getline refline
            getline queryline
            upos = match(refline, /[acgtn]/) - 1
            mid1 = match(refline, /[acgtn]-*[acgtn]/)
            qstart1 = upos + match(substr(queryline, upos + 1), /[ACGTN]/) - 1
            qend1 = qstart1 + match(substr(queryline, qstart1 + 1, mid1 - qstart1), /[ACGTN]-*$/)
            mid2 = mid1 + match(substr(refline, mid1 + 1), /[acgtn]/) - 1
            qstart2 = mid2 + match(substr(queryline, mid2 + 1), /[ACGTN]/) - 1

            seg = substr(refline, 1, qend1)
            rpos1 = gsub(/[acgtnACGTN]/, "", seg)

            seg = substr(queryline, 1, qend1)
            qpos1 = gsub(/[ACGTN]/, "", seg)

            seg = substr(refline, mid2 + 1, qstart2 - mid2)
            rpos2 = gsub(/[acgtnACGTN]/, "", seg)

            seg = substr(queryline, 1, mid2)
            qpos2 = gsub(/[ACGTN]/, "", seg)

            printf("%d\t%d\t%d\t%d\n", rpos1, rpos2, qpos1, qpos2)
        }
    '
}

reverse_fastq() {
    perl -pE '
        if ($. % 4 == 2) {
            tr/acgtACGT/TGCATGCA/;
        }
        if ($. % 2 == 0) {
            chomp;
            $_ = scalar reverse;
            $_ .= "\n";
        }
    '
}

query_to_seq() {
    sq_file=$1
    sort -k1,1 |
    join -t $'\t' -1 2 -2 1 \
        <(sort -k2,2 $sq_file) - |
    cut -f2- |
    sort -k1,1
}

opt_indel() {
    join -t $'\t' -1 1 -2 1 \
        <(
            seq -f "seq%g" $query_per_ref |
            paste - bench/truth.tsv |
            sort -k1,1
        ) - |
    sort -k1,1V |
    utils/compare_indel.py \
        $(cut -f2 bench/ref.ref) \
        $(cut -f5 bench/ref.ref)
}

mkdir -p bench

ref_num=1
query_per_ref=100

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" mode reflen probability software realtime usertime systime memory query rpos1 rpos2 qpos1 qpos2 >bench.tsv
# for mode in single double
# do
#     for ref_len in 100 200 300 400 500
#     do
#         for probability in 0.01 0.02 0.03 0.04 0.05
#         do
for mode in single
do
    for ref_len in 100
    do
        for probability in 0.01
        do
            title "generate random queries"
            utils/random_seq_methods.py $ref_len $ref_num $query_per_ref $probability $mode \
                > bench/ref.ref \
                3> bench/query.post \
                4> bench/truth.tsv
            seq $query_per_ref |
            sed 's/^/seq/' |
            paste - bench/truth.tsv |
            row_print simulate "*" "*" "*" "*" \
                >> bench.tsv

            title "rearr"
            mkdir -p bench/rearr
            # time is a shell keyword, explicitly specify /bin/time to use the time executable
            /bin/time -f '%e\t%U\t%S\t%M' -o bench/rearr/time.linux \
                ../core/Rearrangement/build/rearrangement \
                    < bench/query.post \
                    3< bench/ref.ref \
                    > bench/rearr/query.alg \
                    2> /dev/null
            paste \
                <(seq -f "seq%g" $query_per_ref) \
                bench/truth.tsv \
                <(parse_rqpos < bench/rearr/query.alg) |
            utils/compare_indel.py \
                $(cut -f2 bench/ref.ref) \
                $(cut -f5 bench/ref.ref) |
            row_print rearr $(cat bench/rearr/time.linux) \
                >> bench.tsv

            # title "RESSO"
            # mkdir -p bench/RESSO
            # to_fastq \
            #     < bench/query.post \
            #     > bench/RESSO/query.fq
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/RESSO/time.linux \
            #     CRISPResso \
            #         -r1 bench/RESSO/query.fq \
            #         -a $(get_ref < bench/ref.ref) \
            #         -g $(get_sgRNA < bench/ref.ref) \
            #         -o bench/RESSO \
            #         -amas 0
            # unzip -q -o \
            #     $(find bench/RESSO -name "*.zip") \
            #     -d bench/RESSO
            # seq -f "seq%g" $query_per_ref |
            # paste - <(cut -f1 bench/query.post) \
            #     > bench/RESSO/sq_file
            # utils/ARRANGE_RESULTS.py \
            #     < bench/RESSO/Alleles_frequency_table.txt \
            #     CRISPResso \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) |
            #     query_to_seq bench/RESSO/sq_file |
            #     opt_indel |
            #     row_print RESSO $(cat bench/RESSO/time.linux) \
            #         >> bench.tsv

            # title "AMP"
            # mkdir -p bench/AMP
            # amplicon=$(get_amplicon < bench/ref.ref)
            # primer=$(head -c6 <<<$amplicon | tr 'acgtn' 'ACGTN')
            # to_fastq \
            #     < bench/query.post |
            # sed "2~4{s/^/$primer/}; 4~4{s/^/~~~~~~/}" \
            #     > bench/AMP/query.fq
            # printf "ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor\n" \
            #     > bench/AMP/config.csv
            # printf "ID_1,barcode_1,query.fq,,group_1,0,%s,%s,,0,%s," \
            #     $(get_sgRNA < bench/ref.ref) \
            #     $primer \
            #     $amplicon \
            #     >> bench/AMP/config.csv
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/AMP/time.linux \
            #     utils/AMP.r \
            #         bench/AMP/config.csv \
            #         bench/AMP/ \
            #         bench/AMP/
            # sed -nr '1~4s/^@//; N; s/\n/\t/p' \
            #     < bench/AMP/query.fq \
            #     > bench/AMP/sq_file
            # utils/ARRANGE_RESULTS.py \
            #     < bench/AMP/alignments/alignments.txt \
            #     amplican \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) \
            #     $(expr $(cut -f1 bench/ref.ref) - 6) \
            #     ${#primer} |
            #     paste \
            #         <(
            #             sed -n '2~4p' bench/AMP/query.fq |
            #             sort |
            #             uniq -c |
            #             sort -s -k1,1nr |
            #             awk '{print $2}'
            #         ) - |
            #     query_to_seq bench/AMP/sq_file |
            #     opt_indel |
            #     row_print AMP $(cat bench/AMP/time.linux) \
            #         >> bench.tsv

            # title "CRVS"
            # mkdir -p bench/CRVS
            # get_ref \
            #     < bench/ref.ref |
            # sed '1i >ref' \
            #     > bench/CRVS/ref.fa
            # to_fastq \
            #     < bench/query.post \
            #     > bench/CRVS/query.fq
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/CRVS/time.linux \
            #     bash -c "
            #         bwa index -p bench/CRVS/ref bench/CRVS/ref.fa
            #         bwa mem -v 3 -T -9999 bench/CRVS/ref bench/CRVS/query.fq |
            #         samtools sort -o bench/CRVS/query.bam
            #         samtools index -b bench/CRVS/query.bam
            #         utils/variants.r \
            #             $(sed '1d' bench/CRVS/ref.fa) \
            #             bench/CRVS/query.bam \
            #             $(expr $(cut -f3 bench/ref.ref) - 17) 20 \
            #             > bench/CRVS/result
            #     "
            # cut -f1 bench/CRVS/result |
            # sort |
            # join -t $'\t' -1 1 -2 1 - \
            #     <(
            #         samtools view bench/CRVS/query.bam |
            #         sort -k1,1
            #     ) |
            # utils/ARRANGE_RESULTS.py \
            #     CrisprVariants \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) |
            # opt_indel |
            # row_print CRVS $(cat bench/CRVS/time.linux) \
            #     >> bench.tsv

            # title "ADIV"
            # mkdir -p bench/ADIV
            # to_fastq \
            #     < bench/query.post |
            # sed '2~4{s/^/AAAAAATGTAAAACGACGGCCAGT/; s/$/aagacac/}; 4~4{s/^/~~~~~~~~~~~~~~~~~~~~~~~~/; s/$/~~~~~~~/}' \
            #     > bench/ADIV/query.for.fq
            # reverse_fastq \
            #     < bench/ADIV/query.for.fq \
            #     > bench/ADIV/query.rev.fq
            # printf "plate0\tA0\trefin_ref_45-200\tM\t0\tAAAAAA\n" \
            #     > bench/ADIV/query.bar
            # get_ref \
            #     < bench/ref.ref |
            # sed '1i >ref' \
            #     > bench/ADIV/ref.fa
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/ADIV/time.linux \
            #     bash -c "
            #         novoindex bench/ADIV/ref.nix bench/ADIV/ref.fa
            #         novoalign -t 0,6 -o SAM -d bench/ADIV/ref.nix -f bench/ADIV/query.for.fq bench/ADIV/query.rev.fq |
            #         samtools view -o bench/ADIV/query.bam
            #         rm -rf bench/ADIV/Workdir_target_ bench/ADIV/Output_target_
            #         cd bench/ADIV
            #         ../../tools/ampliconDIVider/ampliconDIVider_driver.sh -b query.bar -f ref.fa -x ref.nix query.bam
            #         cd ../..
            #     "
            # samtools view bench/ADIV/Workdir_target_/query.aligned.bam |
            # utils/ARRANGE_RESULTS.py \
            #     "ampliconDIVider" \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) |
            # sort -k1,1 |
            # opt_indel |
            # row_print ADIV $(cat bench/ADIV/time.linux) \
            #     >> bench.tsv

            # title "CRGR"
            # mkdir -p bench/CRGR
            # to_fastq \
            #     < bench/query.post \
            #     > bench/CRGR/query.for.fq
            # reverse_fastq \
            #     < bench/CRGR/query.for.fq \
            #     > bench/CRGR/query.rev.fq
            # get_ref \
            #     < bench/ref.ref |
            # sed '1i >ref' \
            #     > bench/CRGR/ref.fa
            # rm -rf bench/CRGR/output
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/CRGR/time.linux \
            #     tools/CRISPR-GRANT/bin/indel_analysis \
            #         -1:bench/CRGR/query.for.fq \
            #         -2:bench/CRGR/query.rev.fq \
            #         -r:bench/CRGR/ref.fa \
            #         -o:bench/CRGR/output
            # samtools sort -n bench/CRGR/output/4.Mapping_sorted.bam |
            # samtools view |
            # utils/ARRANGE_RESULTS.py \
            #     "CRISPR-GRANT" \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) |
            # sort -k1,1 |
            # opt_indel |
            # row_print CRGR $(cat bench/CRGR/time.linux) \
            #     >> bench.tsv

            # title "ZhangFeng"
            # mkdir -p bench/ZhangFeng
            # to_fastq \
            #     < bench/query.post \
            #     > bench/ZhangFeng/query.fq
            # printf "query,bench/ZhangFeng/query.fq,%s,%s\n" \
            #     $(get_sgRNA < bench/ref.ref) \
            #     $(get_ref_ZhangFeng < bench/ref.ref) \
            #     > bench/ZhangFeng/samplefile
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/ZhangFeng/time.linux \
            #     tools/Screening_Protocols_manuscript/my_calculate_indel.py \
            #         -i bench/ZhangFeng/samplefile \
            #         -no-m \
            #         -o bench/ZhangFeng/calc_indel_out.csv \
            #         > bench/ZhangFeng/result
            # tail -n+2 bench/ZhangFeng/result |
            # utils/ARRANGE_RESULTS.py \
            #     ZhangFeng \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) \
            #     $(cut -f1 bench/ref.ref) |
            # sort -k1,1 |
            # opt_indel |
            # row_print ZhangFeng $(cat bench/ZhangFeng/time.linux) \
            #     >> bench.tsv

            # title "SelfTarget"
            # mkdir -p bench/SelfTarget
            # printf ">ref %d FORWARD\n%s\n" \
            #     $(expr $(cut -f3 bench/ref.ref) + 3) \
            #     $(get_ref < bench/ref.ref) \
            #     > bench/SelfTarget/oligo.fa
            # to_fastq \
            #     < bench/query.post \
            #     > bench/SelfTarget/query.fq
            # /bin/time -f '%e\t%U\t%S\t%M' -o bench/SelfTarget/time.linux \
            #     tools/SelfTarget-master/indel_analysis/indelmap/indelmap \
            #     bench/SelfTarget/query.fq \
            #     bench/SelfTarget/oligo.fa \
            #     bench/SelfTarget/outputfile \
            #     0
            # tail -n+2 bench/SelfTarget/outputfile |
            # utils/ARRANGE_RESULTS.py \
            #     SelfTarget \
            #     $(cut -f3 bench/ref.ref) \
            #     $(cut -f4 bench/ref.ref) |
            # sort -k1,1 |
            # opt_indel |
            # row_print SelfTarget $(cat bench/SelfTarget/time.linux) \
            #     >> bench.tsv

            title "Indel_searcher_2"
            mkdir -p bench/Indel_searcher_2
            user=ljw
            project=rearr
            cp tools/CRISPR_toolkit/EDNAFULL bench
            cp tools/CRISPR_toolkit/Indel_searcher_2/my_Run_indel_searcher.py bench/Indel_searcher_2
            cp tools/CRISPR_toolkit/Indel_searcher_2/my_Indel_searcher_crispresso_hash.py bench/Indel_searcher_2
            mkdir -p bench/Indel_searcher_2/Core
            cp -r tools/CRISPR_toolkit/Core/__init__.py bench/Indel_searcher_2/Core
            cp -r tools/CRISPR_toolkit/Core/my_CoreSystem.py bench/Indel_searcher_2/Core
            mkdir -p bench/Indel_searcher_2/Input/${user}/FASTQ/${project}/query
            barcode=$(
                cut -f2 bench/ref.ref |
                cut -c1-19
            )
            sed 's/^/'$barcode'/' \
                < bench/query.post |
            to_fastq \
                > bench/Indel_searcher_2/Input/${user}/FASTQ/${project}/query/query.fastq
            mkdir -p bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref
            echo $barcode \
                > bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref/Barcode.txt
            get_ref \
                < bench/ref.ref \
                > bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref/Reference_sequence.txt
            get_ref \
                < bench/ref.ref |
            cut -c20- \
                > bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref/Target_region.txt
            mkdir -p bench/Indel_searcher_2/User/${user}
            printf "query\tref\n" \
                > bench/Indel_searcher_2/User/${user}/${project}.txt
            /bin/time -f '%e\t%U\t%S\t%M' -o bench/Indel_searcher_2/time.linux \
                bench/Indel_searcher_2/my_Run_indel_searcher.py \
                    --python $(which python) \
                    --user ${user} \
                    --project ${project} \
                    --pam_type Cas9 \
                    --pam_pos Forward \
                    -t 15
            sort -k1,1 \
                < bench/Indel_searcher_2/result |
            utils/ARRANGE_RESULTS.py \
                "Indel_searcher_2" \
                $(cut -f3 bench/ref.ref) \
                $(cut -f4 bench/ref.ref) \
                19 |
            opt_indel |
            row_print Indel_searcher_2 $(cat bench/Indel_searcher_2/time.linux) \
                >> bench.tsv
        done
    done
done