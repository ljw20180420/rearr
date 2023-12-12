#!/usr/bin/env python
import Bio.Seq, Bio.Align, numpy

def infer_cut(ref, sgRNA):
    # infer the cut point of reference from sgRNA
    pos = ref.find(sgRNA)
    if pos == -1:
        sgRNA = Bio.Seq.Seq(sgRNA).reverse_complement().__str__()
        pos = ref.find(sgRNA)
        if pos == -1:
            raise Exception("cannot locate sgRNA in reference")
    if ref[pos+len(sgRNA)+1:pos+len(sgRNA)+3] == "GG":
        return pos+len(sgRNA)-3, "NGG"
    if ref[pos-3:pos-1] == "CC":
        return pos+3, "CCN"
    if ref[pos+len(sgRNA)-2:pos+len(sgRNA)] == "GG":
        return pos+len(sgRNA)-6, "NGG"
    if ref[pos:pos+2] == "CC":
        return pos+6, "CCN"
    PAMpos = sgRNA.rfind("GG")
    if PAMpos != -1:
        return PAMpos-4, "NGG"
    PAMpos = sgRNA.find("CC")
    if PAMpos != -1:
        return PAMpos+6, "CCN"
    PAMpos = ref.find("GG", pos+len(sgRNA))
    if PAMpos != -1:
        return PAMpos-4, "NGG"
    PAMpos = ref.rfind("CC", 0, pos)
    if PAMpos != -1:
        return PAMpos+6, "CCN"

def try_reverse_complement(countfile, ref1, cut1, NGGCCNtype1, ref2, cut2, NGGCCNtype2, num=10):
    # use the first num reads in countfile to determine whether reference should be reversed and complemented
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1
    with open(countfile, "r") as cf:
        f1f2 = ref1[:cut1] + ref2[cut2:]
        f1r2 = ref1[:cut1] + Bio.Seq.Seq(ref2[:cut2]).reverse_complement().__str__()
        r1f2 = Bio.Seq.Seq(ref1[cut1:]).reverse_complement().__str__() + ref2[cut2:]
        r1r2 = Bio.Seq.Seq(ref1[cut1:]).reverse_complement().__str__() + Bio.Seq.Seq(ref2[:cut2]).reverse_complement().__str__()
        f1f2scores, f1r2scores, r1f2scores, r1r2scores = [], [], [], []
        for _ in range(num):
            line = cf.readline()
            if not line:
                break
            read = line.split("\t")[0]
            f1f2scores.append(aligner.align(read, f1f2).score)
            f1r2scores.append(aligner.align(read, f1r2).score)
            r1f2scores.append(aligner.align(read, r1f2).score)
            r1r2scores.append(aligner.align(read, r1r2).score)
    f1f2score, f1r2score, r1f2score, r1r2score = numpy.mean(f1f2scores), numpy.mean(f1r2scores), numpy.mean(r1f2scores), numpy.mean(r1r2scores)
    idx = numpy.argmax([f1f2score, f1r2score, r1f2score, r1r2score])
    if idx == 2 or idx == 3:
        ref1 = Bio.Seq.Seq(ref1).reverse_complement().__str__()
        cut1 = len(ref1) - cut1
        NGGCCNtype1 = "CCN" if NGGCCNtype1 == "NGG" else "NGG"
    if idx == 1 or idx == 3:
        ref2 = Bio.Seq.Seq(ref2).reverse_complement().__str__()
        cut2 = len(ref2) - cut2
        NGGCCNtype2 = "CCN" if NGGCCNtype2 == "NGG" else "NGG"
    return ref1, cut1, NGGCCNtype1, ref2, cut2, NGGCCNtype2

def get_ref_file(ref1, cut1, ref2, cut2, ref_file):
    # construct reference file
    with open(ref_file, "w") as rf:
        _ = rf.write(f"0\n{ref1}\n{cut1}\n{cut2}\n{ref2}\n{len(ref2)}\n")