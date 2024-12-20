#!/usr/bin/env -S gawk -f

# Usage: correct_micro_homology.AWK -- refFile correctFile
# correct the microhomology to fit the 5' overhang

BEGIN{
    FS = "\t"
    refFile = ARGV[1]
    correctFile = ARGV[2]
    ref_id = 0
    while (getline ref < refFile) {
        getline correct < correctFile
        n = split(ref, ref_arr, "\t")
        split(correct, correct_type[ref_id], "\t")
        ref_accum_len = 0
        for (i = 3; i < n; i += 3) {
            cut1s[ref_id, i / 3] = ref_arr[i]
            cut1s_accum[ref_id, i / 3] = ref_accum_len + cut1s[ref_id, i / 3]
            ref_accum_len += length(ref_arr[i - 1])
            cut2s[ref_id, i / 3] = ref_arr[i + 1]
            cut2s_accum[ref_id, i / 3] = ref_accum_len + cut2s[ref_id, i / 3]
        }
        ++ref_id
    }
    for (i = 1; i <= 2; ++i) {
        delete ARGV[i]
    }
}

function query_patsplit(refs, dashs, queryline, targets, inserts,       start, i) {
    start = 1
    for (i = 1; i <= length(refs); ++i) {
        inserts[i - 1] = substr(queryline, start, length(dashs[i - 1]))
        start += length(dashs[i - 1])
        targets[i] = substr(queryline, start, length(refs[i]))
        start += length(refs[i])
    }
    inserts[length(refs)] = substr(queryline, start, length(dashs[length(refs)]))
}

function query_seg_range(target, seg_range,       mat) {
    seg_range[1] = match(target, /([ACGTN][-ACGTN]*[ACGTN]|[ACGTN])/, mat) - 1
    seg_range[2] = seg_range[1] + length(mat[0])
}

function get_gap_cut(ref, cut,       n, segs, gaps, accum_len, j, k) {
    n = patsplit(ref, segs, /[acgtnACGTN]+/, gaps)
    accum_len = 0
    for (j = 1; j <= n; ++j) {
        accum_len += length(segs[j])
        if (accum_len >= cut) {
            break
        }
    }
    gap_cut = cut
    for (k = 1; k < j; ++k) {
        gap_cut += length(gaps[k])
    }
    return gap_cut
}

function longest_common_fix(string1, string2, fix,      i, start, rgx) {
    string1 = toupper(string1)
    string2 = toupper(string2)
    for (i = 1; i <= length(string2); ++i) {
        start = fix == "prefix" ? 1 : length(string2) - i + 1
        rgx = fix == "prefix" ? "^" substr(string2, start, i) : substr(string2, start, i) "$"
        if (string1 !~ rgx) {
            return i - 1
        }
    }
}

function print_mark(insert, seg_range, ref, target,       ref_block, query_block) {
    printf("%s\t", insert)
    query_pos += length(insert)
    ref_pos += (seg_range[1] == -1 ? 0 : seg_range[1])
    printf("%d\t%d\t", ref_pos, query_pos)
    ref_block = substr(ref, seg_range[1] + 1, seg_range[2] - seg_range[1])
    gsub(/-/, "", ref_block)
    ref_pos += length(ref_block)
    query_block = substr(target, seg_range[1] + 1, seg_range[2] - seg_range[1])
    gsub(/-/, "", query_block)
    query_pos += length(query_block)
    printf("%d\t%d\t", ref_pos, query_pos)
    ref_pos += seg_range[1] == -1 ? length(ref) : length(ref) - seg_range[2]
}

{
    idx = $1
    count = $2
    score = $3
    ref_id = $4
    printf("%d\t%d\t%d\t%d\t", idx, count, score, ref_id)

    getline refline
    getline queryline

    patsplit(refline, refs, /[acgtn][-ACGTN]*[acgtn]/, dashs)
    query_patsplit(refs, dashs, queryline, targets, inserts)

    ref_pos = query_pos = 0
    query_seg_range(targets[1], seg_range1)
    for (i = 1; i < length(refs); ++i) {
        query_seg_range(targets[i + 1], seg_range2)
        if (length(inserts[i]) == 0 && seg_range1[1] != 0 && seg_range2[1] != 0) {
            gap_cut1 = get_gap_cut(refs[i], cut1s[ref_id, i])
            gap_cut2 = get_gap_cut(refs[i + 1], cut2s[ref_id, i])
            
            correct_direct = ""
            if (correct_type[ref_id][i] == "up") {
                if (seg_range1[2] < gap_cut1) {
                    receiver = substr(refs[i], seg_range1[2] + 1, gap_cut1 - seg_range1[2])
                    provider = substr(refs[i + 1], seg_range2[1] + 1, seg_range2[2] - seg_range2[1])
                    correct_direct = "prefix"
                } else if (seg_range1[2] > gap_cut1) {
                    provider = substr(refs[i], gap_cut1 + 1, seg_range1[2] - gap_cut1)
                    receiver = substr(refs[i + 1], 1, seg_range2[1])
                    correct_direct = "suffix"
                }
            } else {
                if (seg_range2[1] < gap_cut2) {
                    provider = substr(refs[i + 1], seg_range2[1] + 1, gap_cut2 - seg_range2[1])
                    receiver = substr(refs[i], seg_range1[2] + 1)
                    correct_direct = "prefix"
                } else if (seg_range2[1] > gap_cut2) {
                    receiver = substr(refs[i + 1], gap_cut2 + 1, seg_range2[1] - gap_cut2)
                    provider = substr(refs[i], seg_range1[1] + 1, seg_range1[2] - seg_range1[1])
                    correct_direct = "suffix"
                }
            }

            if (correct_direct != "") {
                correct_length = longest_common_fix(receiver, provider, correct_direct)
                target1_split_start = correct_direct == "prefix" ? seg_range1[2] : seg_range1[2] - correct_length
                target2_split_start = correct_direct == "prefix" ? seg_range2[1] : seg_range2[1] - correct_length
                target1_exchange = substr(targets[i], target1_split_start + 1, correct_length)
                target2_exchange = substr(targets[i + 1], target2_split_start + 1, correct_length)
                targets[i] = substr(targets[i], 1, target1_split_start) target2_exchange substr(targets[i], target1_split_start + correct_length + 1)
                targets[i + 1] = substr(targets[i + 1], 1, target2_split_start) target1_exchange substr(targets[i + 1], target2_split_start + correct_length + 1)

                shift = correct_direct == "prefix" ? correct_length : -correct_length
                seg_range1[2] += shift
                seg_range2[1] += shift
            }
        }
        print_mark(inserts[i - 1], seg_range1, refs[i], targets[i])
        seg_range1[1] = seg_range2[1]
        seg_range1[2] = seg_range2[2]
    }
    print_mark(inserts[i - 1], seg_range1, refs[i], targets[i])
    printf("%s", inserts[i])
    for (i = 1; i < length(refs); ++i) {
        printf("\t%d\t%d", cut1s_accum[ref_id, i], cut2s_accum[ref_id, i])
    }
    queryline = inserts[0]
    for (i = 1; i <= length(refs); ++i) {
        queryline = queryline targets[i] inserts[i]
    }
    printf("\n%s\n%s\n", refline, queryline)
}
