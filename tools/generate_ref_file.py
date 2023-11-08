#!/usr/bin/env python
import sys, pre_handle_sgRNA_PAM_ref

input, ref, sgRNA = sys.argv[1], sys.argv[2], sys.argv[3]
ext1 = int(sys.argv[4]) if len(sys.argv) > 4 else 30
ext2 = int(sys.argv[5]) if len(sys.argv) > 5 else ext1

ref = pre_handle_sgRNA_PAM_ref.try_reverse_complement(f"{input}.count", ref, num=10) # use the first num reads in countfile to determine whether reference should be reversed and complemented
# ref = pre_handle_sgRNA_PAM_ref.auto_extend_ref(f"{input}.count", ref, num=100) # use the first num reads in countfile to determine how long ref should be extended
cut, NGGCCNtype = pre_handle_sgRNA_PAM_ref.infer_cut(ref, sgRNA) # infer the cut point of reference from sgRNA
pre_handle_sgRNA_PAM_ref.get_ref_file(ref, cut, f"{input}.ref.{cut}.{ext1}.{ext2}", ext1=ext1, ext2=ext2) # construct reference file
sys.stdout.write(f"{cut} {NGGCCNtype}")