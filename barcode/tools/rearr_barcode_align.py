#!/usr/bin/env python
import sys, os, subprocess
_, fq1, ext1up, ext2up = sys.argv
ext1up, ext2up = int(ext1up), int(ext2up)
_ = subprocess.run(f'''> {fq1}.alg''', shell=True, check=True)
_ = subprocess.run(f'''printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "barcode" "index" "count" "score" "updangle" "ref_start1" "query_start1" "ref_end1" "query_end1" "random_insertion" "ref_start2" "query_start2" "ref_end2" "query_end2" "downdangle" "cut1" "cut2" > {fq1}.table''', shell=True, check=True, executable="/bin/bash")
barcodegroup = ""
for bcdline in sys.stdin:
    _, count, naseq, barcode_end, barcode, ref1, ref2 = bcdline.rstrip().split('\t')
    barcode_end = int(barcode_end)
    if barcodegroup != barcode:
        if barcodegroup:
            cf.close()
            _ = subprocess.run(f'''rearrangement <{fq1}.countfile 3<{fq1}.reference -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | correct_micro_homology.py {ext1up} NGG {ext2up} NGG {len(ref1)} | tee -a {fq1}.alg | awk -v OFS="\t" -v barcode={barcodegroup} -v cut1={ext1up} -v cut2={len(ref1)+ext2up} 'NR%3==1{{print barcode, $0, cut1, cut2}}' >> {fq1}.table; echo >> {fq1}.alg''', shell=True, check=True)
        barcodegroup = barcode
        with open(f"{fq1}.reference", "w") as rd:
            _ = rd.write(f"0\n{ref1}\n{ext1up}\n{ext2up}\n{ref2}\n{len(ref2)}\n")
        cf = open(f"{fq1}.countfile", "w")
        _ = subprocess.run(f'''echo "{barcodegroup}" >> {fq1}.alg''', shell=True, check=True)
    _ = cf.write(f"{naseq[barcode_end + 3:]}\t{count}\n")
        
cf.close()
_ = subprocess.run(f'''rearrangement <{fq1}.countfile 3<{fq1}.reference -u -3 -v -9 -s0 -6 -s1 4 -s2 2 -qv -9 | correct_micro_homology.py {ext1up} NGG {ext2up} NGG {len(ref1)} | tee -a {fq1}.alg | awk -v OFS="\t" -v barcode={barcodegroup} -v cut1={ext1up} -v cut2={len(ref1)+ext2up} 'NR%3==1{{print barcode, $0, cut1, cut2}}' >> {fq1}.table; echo >> {fq1}.alg''', shell=True, check=True)
os.remove(f"{fq1}.countfile")
os.remove(f"{fq1}.reference")
