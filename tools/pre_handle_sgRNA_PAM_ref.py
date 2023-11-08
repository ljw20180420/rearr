#!/usr/bin/env python
import subprocess, Bio.Seq, Bio.Align, numpy, os, re, sys, correct_micro_homology

def extend_ref(ref, upstream_extend=0, downstream_extend=0):
    # extend ref to up/downstream if ref comes from genome
    if not upstream_extend and not downstream_extend:
        return ref
    _, flag, chr, pos, left = subprocess.run(f"echo {ref} | bowtie2 -x /home/ljw/hg19_with_bowtie2_index/hg19 -r -U - | samtools view", shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().split('\t', 4)
    mats = re.search("AS:i:(-?\d+)", left)
    if not mats:
        sys.stderr.write("the reference cannot be found in hg19\n")
        return ref
    if int(mats.group(1)) != 0:
        sys.stderr.write("warning: the reference is not exactly match the genome\n")
    pos = int(pos) - 1
    strand = "-" if (int(flag) / 16) % 2 else "+"
    if strand == "-":
        startup, endup, startdown, enddown = pos+len(ref), pos+len(ref)+upstream_extend, pos-downstream_extend, pos
    else:
        startup, endup, startdown, enddown = pos-upstream_extend, pos, pos+len(ref), pos+len(ref)+downstream_extend
    extup = subprocess.run(f''' printf "{chr}\t{startup}\t{endup}\t.\t.\t{strand}\n" | bedtools getfasta -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed - | tail -n1 ''', shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().rstrip()
    extdown = subprocess.run(f''' printf "{chr}\t{startdown}\t{enddown}\t.\t.\t{strand}\n" | bedtools getfasta -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed - | tail -n1 ''', shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().rstrip()
    return (extup + ref + extdown).upper()

def auto_extend_ref(countfile, ref, num=100):
    # use the first num reads in countfile to determine how long ref should be extended
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1
    ref = extend_ref(ref, upstream_extend=100, downstream_extend=100)
    with open(countfile, "r") as cf:
        start, end = len(ref), 0
        for _ in range(num):
            line = cf.readline()
            if not line:
                break
            read = line.split("\t")[0]
            ranges = aligner.align(read, ref)[0].coordinates[1]
            start = min(start, ranges[0])
            end = max(end, ranges[-1])
    start = max(0, start-10)
    end = min(len(ref), end+10)
    return ref[start:end]

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

def try_reverse_complement(countfile, ref, num=10):
    # use the first num reads in countfile to determine whether reference should be reversed and complemented
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1
    with open(countfile, "r") as cf:
        revref = Bio.Seq.Seq(ref).reverse_complement().__str__()
        scores, revscores = [], []
        for _ in range(num):
            line = cf.readline()
            if not line:
                break
            read = line.split("\t")[0]
            scores.append(aligner.align(read, ref).score)
            revscores.append(aligner.align(read, revref).score)
    if numpy.mean(revscores) > numpy.mean(scores):
        return revref
    return ref

def get_ref_file(ref, cut, ref_file, ext1=10, ext2=10):
    # construct reference file
    with open(ref_file, "w") as rf:
        _ = rf.write(f"0\t0\t0\n{ref[:cut + ext1]}\n{cut}\t{cut}\t{cut}\n{ext2}\t{ext2}\t{ext2}\n{ref[cut-ext2:]}\n{len(ref)-cut+ext2}\t{len(ref)-cut+ext2}\t{len(ref)-cut+ext2}\n")

# this method is problematic because extending too long at the cut point leads to abnormal alignments, so do not use it
def auto_cut_extend(countfile, ref, cut, exec, ref_file, num=100):
    # use the first num reads in countfile to determine how long the cut pos should be extended to hold the templated insertion
    ext1 = min(100, len(ref) - cut)
    ext2 = min(100, cut)
    get_ref_file(ref, cut, "tempref", ext1=ext1, ext2=ext2)
    output = subprocess.run(f'''{exec} -file <(head -n{num} {countfile}) -ref_file tempref -ALIGN_MAX 1 -THR_NUM 1 -u -1 -v -3 -s0 -2 -qv -3''', shell=True,  executable="/bin/bash", check=True, capture_output=True).stdout.decode().splitlines()
    os.remove("tempref")
    ext1correct, ext2correct = 0, 0
    for header, _, _ in correct_micro_homology.correct_micro(cut + ext1, cut, cut + ext1 + ext2, output):
        _, _, _, _, ue, _, ds, _ = header.split("\t")
        ext1correct = max(ext1correct, int(ue) - cut)
        ext2correct = max(ext2correct, cut + ext1 + ext2 - int(ds))
    ext1correct = min(ext1correct+5, ext1)
    ext2correct = min(ext2correct+5, ext2)
    get_ref_file(ref, cut, ref_file, ext1=ext1correct, ext2=ext2correct)
    return ext1correct, ext2correct