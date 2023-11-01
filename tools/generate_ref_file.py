#!/usr/bin/env python
import sys, pre_handle_sgRNA_PAM_ref

countfile, ref, sgRNA, exec = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
ext1 = int(sys.argv[5]) if len(sys.argv) > 5 else 10
ext2 = int(sys.argv[6]) if len(sys.argv) > 6 else ext1

ref = pre_handle_sgRNA_PAM_ref.try_reverse_complement(countfile, ref, num=10)
ref = pre_handle_sgRNA_PAM_ref.auto_extend_ref(countfile, ref, num=100)
cut = pre_handle_sgRNA_PAM_ref.infer_cut(ref, sgRNA)
pre_handle_sgRNA_PAM_ref.get_ref_file(ref, cut, "ref_file", ext1=ext1, ext2=ext2)
sys.stdout.write(f"{cut}")