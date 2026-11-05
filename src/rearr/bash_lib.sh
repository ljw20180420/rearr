bench_rearr() {
	local output=$1
	local timetab=$2
	local ref_file=$3
	local query_file=$4
	local alg_file=$5

	# time is a shell keyword, use /bin/time to use the binary
	/bin/time -f '%e\t%U\t%S\t%M' \
		-o "${output}" \
		rearrangement \
		<"${query_file}" \
		3<"${ref_file}" \
		>"${alg_file}" \
		2>"${timetab}"
}
