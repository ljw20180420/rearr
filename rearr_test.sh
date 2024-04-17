#!/bin/bash
rearrangement <test/random.fq.count 3<test/random.fq.ref -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | gawk -f correct_micro_homology.AWK -- test/random.fq.ref NGG NGG >test/random.fq.alg