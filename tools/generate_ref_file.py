#!/usr/bin/env python
import sys, pre_handle_sgRNA_PAM_ref

input, ref1, ref2, sgRNA1, sgRNA2 = sys.argv[1:6]

cut1, NGGCCNtype1 = pre_handle_sgRNA_PAM_ref.infer_cut(ref1, sgRNA1) # infer the cut point of reference from sgRNA
cut2, NGGCCNtype2 = pre_handle_sgRNA_PAM_ref.infer_cut(ref2, sgRNA2) # infer the cut point of reference from sgRNA

ref1, cut1, NGGCCNtype1, ref2, cut2, NGGCCNtype2 = pre_handle_sgRNA_PAM_ref.try_reverse_complement(f"{input}.count", ref1, cut1, NGGCCNtype1, ref2, cut2, NGGCCNtype2, num=10) # use the first num reads in countfile to determine whether reference should be reversed and complemented

pre_handle_sgRNA_PAM_ref.get_ref_file(ref1, cut1, ref2, cut2, f"{input}.ref.{cut1}.{cut2}.{len(ref1)}") # construct reference file
sys.stdout.write(f"{cut1} {NGGCCNtype1} {cut2} {NGGCCNtype2}")