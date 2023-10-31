import os, sys, subprocess

def extend_ref(ref, upstream_extend=0, downstream_extend=0):
    if not upstream_extend and not downstream_extend:
        return ref
    _, flag, chr, pos, _ = subprocess.run(f"echo {ref} | bowtie2 -x /home/ljw/hg19_with_bowtie2_index/hg19 -r -U - | samtools view", shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().split('\t', 4)
    pos = int(pos) - 1
    strand = "-" if (int(flag) / 16) % 2 else "+"
    if strand == "-":
        startup, endup, startdown, enddown = pos+len(ref), pos+len(ref)+upstream_extend, pos-downstream_extend, pos
    else:
        startup, endup, startdown, enddown = pos-upstream_extend, pos, pos+len(ref), pos+len(ref)+downstream_extend
    extup = subprocess.run(f''' printf "{chr}\t{startup}\t{endup}\t.\t.\t{strand}\n" | bedtools getfasta -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed - | tail -n1 ''', shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().rstrip()
    extdown = subprocess.run(f''' printf "{chr}\t{startdown}\t{enddown}\t.\t.\t{strand}\n" | bedtools getfasta -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed - | tail -n1 ''', shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().rstrip()
    return extup + ref + extdown

def get_ref_file(ref, cut, ext=0, reverse_complement=False, ref_file="ljhlyz/ref_file"):
    if reverse_complement:
        ref = subprocess.run(f"echo {ref} | rev | tr 'ACGT' 'TGCA'", shell=True, executable="/bin/bash", check=True, capture_output=True).stdout.decode().rstrip()
        cut = len(ref) - cut
    with open(ref_file, "w") as rf:
        _ = rf.write(f"0\t0\t0\n{ref[:cut + ext]}\n{cut}\t{cut}\t{cut}\n{ext}\t{ext}\t{ext}\n{ref[cut-ext:]}\n{len(ref)-cut+ext}\t{len(ref)-cut+ext}\t{len(ref)-cut+ext}\n")
    return cut



fqgz, ref, cut, ext, reverse_complement = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), bool(sys.argv[5])
ref = extend_ref(ref, upstream_extend=30, downstream_extend=0)
ref_file = "ljhlyz/ref_file"
cut = get_ref_file(ref, cut, ext=ext, reverse_complement=reverse_complement, ref_file="ljhlyz/ref_file")

subprocess.run(f'''zcat {fqgz} | sed -n '2~4p' | sort | uniq -c | sort -k1,1nr | awk -v OFS="\t" '{{print $2, $1}}' >{fqgz}.count''', shell=True, executable="/bin/bash", check=True)

_ = subprocess.run(f'''Rearrangement/build/rearrangement -file {fqgz}.count -ref_file {ref_file} -ALIGN_MAX 1 -THR_NUM 24 -u -1 -v -3 -s0 -2 -qv -3 | sed -nr 'N;N;s/\\n/\\t/g;p' | sort -k1,1n | awk -F "\t" '{{for (i=1; i<=NF-3; ++i) printf("%s\\t",$i); printf("%s\\n%s\\n%s\\n", $(NF-2), $(NF-1), $NF);}}' | tee {fqgz}.alg | python barcode/correct_micro_homology.py {cut + ext} {cut} >{fqgz}.correct''', shell=True,  executable="/bin/bash", check=True)

_ = subprocess.run(f'''cat {fqgz}.correct | head -n150 | python ljhlyz/align_align.py {cut + ext} {cut} {cut + ext + ext} >{fqgz}.algalg''', shell=True,  executable="/bin/bash", check=True)