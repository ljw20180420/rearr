#!/bin/bash

# Change directory to the script directory.
cd $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

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
for mode in single; do
	for ref_len in 100; do
		for probability in 0.01; do
			title "ADIV"
			mkdir -p bench/ADIV
			to_fastq \
				<bench/query.post |
				sed '2~4{s/^/AAAAAATGTAAAACGACGGCCAGT/; s/$/aagacac/}; 4~4{s/^/~~~~~~~~~~~~~~~~~~~~~~~~/; s/$/~~~~~~~/}' \
					>bench/ADIV/query.for.fq
			reverse_fastq \
				<bench/ADIV/query.for.fq \
				>bench/ADIV/query.rev.fq
			printf "plate0\tA0\trefin_ref_45-200\tM\t0\tAAAAAA\n" \
				>bench/ADIV/query.bar
			get_ref \
				<bench/ref.ref |
				sed '1i >ref' \
					>bench/ADIV/ref.fa
			/bin/time -f '%e\t%U\t%S\t%M' -o bench/ADIV/time.linux \
				bash -c "
			        novoindex bench/ADIV/ref.nix bench/ADIV/ref.fa
			        novoalign -t 0,6 -o SAM -d bench/ADIV/ref.nix -f bench/ADIV/query.for.fq bench/ADIV/query.rev.fq |
			        samtools view -o bench/ADIV/query.bam
			        rm -rf bench/ADIV/Workdir_target_ bench/ADIV/Output_target_
			        cd bench/ADIV
			        ../../tools/ampliconDIVider/ampliconDIVider_driver.sh -b query.bar -f ref.fa -x ref.nix query.bam
			        cd ../..
			    "
			samtools view bench/ADIV/Workdir_target_/query.aligned.bam |
				utils/ARRANGE_RESULTS.py \
					"ampliconDIVider" \
					$(cut -f3 bench/ref.ref) \
					$(cut -f4 bench/ref.ref) |
				sort -k1,1 |
				opt_indel |
				row_print ADIV $(cat bench/ADIV/time.linux) \
					>>bench.tsv

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

			# title "Indel_searcher_2"
			# mkdir -p bench/Indel_searcher_2
			# user=ljw
			# project=rearr
			# cp tools/CRISPR_toolkit/EDNAFULL bench
			# cp tools/CRISPR_toolkit/Indel_searcher_2/my_Run_indel_searcher.py bench/Indel_searcher_2
			# cp tools/CRISPR_toolkit/Indel_searcher_2/my_Indel_searcher_crispresso_hash.py bench/Indel_searcher_2
			# mkdir -p bench/Indel_searcher_2/Core
			# cp -r tools/CRISPR_toolkit/Core/__init__.py bench/Indel_searcher_2/Core
			# cp -r tools/CRISPR_toolkit/Core/my_CoreSystem.py bench/Indel_searcher_2/Core
			# mkdir -p bench/Indel_searcher_2/Input/${user}/FASTQ/${project}/query
			# barcode=$(
			#     cut -f2 bench/ref.ref |
			#     cut -c1-19
			# )
			# sed 's/^/'$barcode'/' \
			#     < bench/query.post |
			# to_fastq \
			#     > bench/Indel_searcher_2/Input/${user}/FASTQ/${project}/query/query.fastq
			# mkdir -p bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref
			# echo $barcode \
			#     > bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref/Barcode.txt
			# get_ref \
			#     < bench/ref.ref \
			#     > bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref/Reference_sequence.txt
			# get_ref \
			#     < bench/ref.ref |
			# cut -c20- \
			#     > bench/Indel_searcher_2/Input/${user}/Reference/${project}/ref/Target_region.txt
			# mkdir -p bench/Indel_searcher_2/User/${user}
			# printf "query\tref\n" \
			#     > bench/Indel_searcher_2/User/${user}/${project}.txt
			# /bin/time -f '%e\t%U\t%S\t%M' -o bench/Indel_searcher_2/time.linux \
			#     bench/Indel_searcher_2/my_Run_indel_searcher.py \
			#         --python $(which python) \
			#         --user ${user} \
			#         --project ${project} \
			#         --pam_type Cas9 \
			#         --pam_pos Forward \
			#         -t 15
			# sort -k1,1 \
			#     < bench/Indel_searcher_2/result |
			# utils/ARRANGE_RESULTS.py \
			#     "Indel_searcher_2" \
			#     $(cut -f3 bench/ref.ref) \
			#     $(cut -f4 bench/ref.ref) \
			#     19 |
			# opt_indel |
			# row_print Indel_searcher_2 $(cat bench/Indel_searcher_2/time.linux) \
			#     >> bench.tsv

			# title "CRISPR-A"
			# mkdir -p bench/CRISPR-A
			# to_fastq \
			#     < bench/query.post \
			#     > bench/CRISPR-A/query.fq
			# get_ref \
			#     < bench/ref.ref |
			# sed '1i >ref' \
			#     > bench/CRISPR-A/ref.fa
			# /bin/time -f '%e\t%U\t%S\t%M' -o bench/CRISPR-A/time.linux \
			#     bash -c "
			#         tools/crispr-a_nextflow/bin/revComp-fasta.R \
			#             bench/CRISPR-A/ref.fa \
			#             bench/CRISPR-A/ref.correctOrient.fa \
			#             $(get_sgRNA < bench/ref.ref)
			#         minimap2 \
			#             -d bench/CRISPR-A/ref.correctOrient.mmi \
			#             bench/CRISPR-A/ref.correctOrient.fa
			#         minimap2 -a -A 29 -B 17 -O 25 -E 2 \
			#             bench/CRISPR-A/ref.correctOrient.mmi \
			#             bench/CRISPR-A/query.fq \
			#             > bench/CRISPR-A/query.sam
			#     "
			# samtools view \
			#     bench/CRISPR-A/query.sam |
			# sort -k1,1 |
			# utils/ARRANGE_RESULTS.py \
			#     "CRISPR-A" \
			#     $(cut -f3 bench/ref.ref) \
			#     $(cut -f4 bench/ref.ref) |
			# opt_indel |
			# row_print "CRISPR-A" $(cat bench/CRISPR-A/time.linux) \
			#     >> bench.tsv

			# title "SIQ"
			# mkdir -p bench/SIQ
			# to_fastq \
			#     < bench/query.post \
			#     > bench/SIQ/query.fq
			# get_ref \
			#     < bench/ref.ref |
			# sed '1i >ref' \
			#     > bench/SIQ/ref.fa
			# left_flank=$(
			#     cut -f2 bench/ref.ref |
			#     cut -c36-50
			# )
			# right_flank=$(
			#     cut -f5 bench/ref.ref |
			#     cut -c51-65
			# )
			# /bin/time -f '%e\t%U\t%S\t%M' -o bench/SIQ/time.linux \
			#     java my.my_SIQ \
			#         bench/SIQ/query.fq \
			#         bench/SIQ/ref.fa \
			#         ${left_flank} \
			#         ${right_flank} \
			#         bench/SIQ/output.tsv
			# sort -k1,1 bench/SIQ/output.tsv |
			# utils/ARRANGE_RESULTS.py \
			#     "SIQ" \
			#     $(cut -f3 bench/ref.ref) \
			#     $(cut -f4 bench/ref.ref) |
			# opt_indel |
			# row_print "SIQ" $(cat bench/SIQ/time.linux) \
			#     >> bench.tsv
		done
	done
done
