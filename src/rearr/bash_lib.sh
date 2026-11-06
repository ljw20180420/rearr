bench_resso() {
	local linux_file=$1
	local ref=$2
	local sgRNA=$3
	local query_file=$4
	local output_dir=$5

	/bin/time -f '%e,%U,%S,%M' \
		-o "${linux_file}" \
		CRISPResso \
		-r1 "${query_file}" \
		-a "${ref}" \
		-amas 0 \
		--default_min_aln_score 0 \
		-g "${sgRNA}" \
		-o "${output_dir}" \
		--suppress_report \
		--suppress_plots
}
