#!/usr/bin/env python

def get_mmej(ref1_end, ref2_start, ref1, ref2, cut1=50, cut2=10):
    # handle the bug of ref1_end == len(ref1)
    i = ref1_end
    for i in range(ref1_end, cut1 + 1):
        if i >= cut1 or ref2_start + i - ref1_end >= len(ref2):
            break
        if ref1[i] != ref2[ref2_start + i - ref1_end]:
            break
    # handle the bug of ref2_start == 0
    j = ref2_start - 1
    for j in range(ref2_start - 1, cut2 - 2, -1):
        if j < cut2 or ref1_end - ref2_start + j < 0:
            break
        if ref2[j] != ref1[ref1_end - ref2_start + j]:
            break
    # breakpoint()
    return (i - ref1_end) + (ref2_start - j - 1), ref2[j + 1:ref2_start] + ref1[ref1_end:i]

if __name__ == "__main__":
    ref1 = "CTTCACGCGGCGCATGCAGAAGGGCCTGGGGTGGAAGGCGTCACCGTACT"
    ref2 = "TTACCGCACTTCACGGAGATCTTGCGATGCAGGGCGGGGCTCATCTCGTGTGGCAGCTTGGCCATGCTGAAGAGCACGTAGAACATCTCGTCCACGTTGGTGTTCTTCTT"
    ref1_end = 3
    ref2_start = 11
    count, bases = get_mmej(ref1_end, ref2_start, ref1, ref2, cut2=5)
    print(count)
    print(bases)



