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

query_to_seq() {
	sq_file=$1
	sort -k1,1 |
		join -t $'\t' -1 2 -2 1 \
			<(sort -k2,2 $sq_file) - |
		cut -f2- |
		sort -k1,1
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
			title "CRISPR-A"
			mkdir -p bench/CRISPR-A
			to_fastq \
				<bench/query.post \
				>bench/CRISPR-A/query.fq
			get_ref \
				<bench/ref.ref |
				sed '1i >ref' \
					>bench/CRISPR-A/ref.fa
			/bin/time -f '%e\t%U\t%S\t%M' -o bench/CRISPR-A/time.linux \
				bash -c "
			        tools/crispr-a_nextflow/bin/revComp-fasta.R \
			            bench/CRISPR-A/ref.fa \
			            bench/CRISPR-A/ref.correctOrient.fa \
			            $(get_sgRNA <bench/ref.ref)
			        minimap2 \
			            -d bench/CRISPR-A/ref.correctOrient.mmi \
			            bench/CRISPR-A/ref.correctOrient.fa
			        minimap2 -a -A 29 -B 17 -O 25 -E 2 \
			            bench/CRISPR-A/ref.correctOrient.mmi \
			            bench/CRISPR-A/query.fq \
			            > bench/CRISPR-A/query.sam
			    "
			samtools view \
				bench/CRISPR-A/query.sam |
				sort -k1,1 |
				utils/ARRANGE_RESULTS.py \
					"CRISPR-A" \
					$(cut -f3 bench/ref.ref) \
					$(cut -f4 bench/ref.ref) |
				opt_indel |
				row_print "CRISPR-A" $(cat bench/CRISPR-A/time.linux) \
					>>bench.tsv

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
