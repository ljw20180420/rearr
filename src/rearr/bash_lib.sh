bench_rearr() {
	local linux_file=$1
	local time_file=$2
	local ref_file=$3
	local query_file=$4
	local alg_file=$5

	# time is a shell keyword, use /bin/time to use the binary
	/bin/time -f '%e,%U,%S,%M' \
		-o "${linux_file}" \
		core/Rearrangement/build/rearrangement \
		<"${query_file}" \
		3<"${ref_file}" \
		>"${alg_file}" \
		4>"${time_file}"
}

bench_resso() {
	local linux_file=$1
	local ref=$3
	local sgRNA=$4
	local query_file=$5
	local output_dir=$6

	/bin/time -f '%e,%U,%S,%M' \
		-o "${linux_file}" \
		CRISPResso \
		-r1 "${query_file}" \
		-a "${ref}" \
		-g "${sgRNA}" \
		-o "${output_dir}" \
		-amas 0 \
		4>"${time_file}"
}
