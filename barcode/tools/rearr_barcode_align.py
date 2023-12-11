#!/usr/bin/env python
import sys, os, subprocess, Bio.Seq, more_itertools
_, barcodefile, csvfile, bowtie2genome, getfastagenome = sys.argv
# scaffold = Bio.Seq.Seq("gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc").reverse_complement().__str__().upper()
extup_left, extdown_left, extup_right, extdown_right = 50, 0, 10, 100
min_seq = 20
with open(barcodefile, "r") as bcd, os.popen(f'''perl -anF, -E 'if($.>1){{substr($F[9], 16, 2)="CC"; $F[11]=~tr/ACGT/TGCA/; say $F[9] , "\t" , scalar reverse $F[11]}}' {csvfile} | sort -k2,2 | cut -f1 | bowtie2 --quiet -x {bowtie2genome} -r -U - 2> /dev/null | samtools view | rearr_barcode2ref.py {extup_left} {extdown_left} {extup_right} {extdown_right} | bedtools getfasta -s -fi {getfastagenome} -bed - | sed '1~2d' ''', "r") as bf, os.popen(f'''tail -n+2 {csvfile} | cut -d',' -f12 | tr ACGT TGCA | rev | sort ''', "r") as rcb:
    _ = subprocess.run(f'''> {barcodefile}.alg''', shell=True, check=True)
    _ = subprocess.run(f'''printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "barcode" "index" "count" "score" "updangle" "ref_start1" "query_start1" "ref_end1" "query_end1" "random_insertion" "ref_start2" "query_start2" "ref_end2" "query_end2" "downdangle" "cut1" "cut2" > {barcodefile}.table''', shell=True, check=True, executable="/bin/bash")

    bcdline = bcd.readline().rstrip()
    if not bcdline:
        raise Exception("*.fq/fastq.barcode is empty")
    _, count, naseq, barcode_end, barcodegroup = bcdline.split('\t')
    barcode_end = int(barcode_end)
    for ref_left, ref_right in more_itertools.batched(bf, 2):
        ref_left, ref_right = ref_left.rstrip(), ref_right.rstrip()
        barcode = rcb.readline().rstrip()
        if barcode != barcodegroup:
            continue
        cf = open(f"{barcodefile}.countfile", "w")
        _ = cf.write(f"{naseq[barcode_end + 3:]}\t{count}\n")
        _ = subprocess.run(f'''echo "{barcodegroup}" >> {barcodefile}.alg''', shell=True, check=True)
        with open(f"{barcodefile}.reference", "w") as rd:
            _ = rd.write(f"0\n{ref_left}\n{extup_left}\n{extup_right}\n{ref_right}\n{extup_right+extdown_right}\n")
        
        while True:
            bcdline = bcd.readline().rstrip()
            if not bcdline:
                break
            _, count, naseq, barcode_end, barcode = bcdline.split('\t')
            barcode_end = int(barcode_end)
            if barcode != barcodegroup:
                barcodegroup = barcode
                break
            _ = cf.write(f"{naseq[barcode_end + 3:]}\t{count}\n")
        cf.close()
        runres = subprocess.run(f'''rearrangement <{barcodefile}.countfile 3<{barcodefile}.reference -u -1 -v -3 -s0 -2 -qv -3 | correct_micro_homology.py {extup_left} {extdown_left} {extup_right} NGG | tee -a {barcodefile}.alg | awk -v OFS="\t" -v barcode={barcodegroup} -v cut1={extup_left} -v cut2={len(ref_left)+extup_right} 'NR%3==1{{print barcode, $0, cut1, cut2}}' >> {barcodefile}.table; echo >> {barcodefile}.alg''', shell=True, check=True)
    cf.close()
    runres = subprocess.run(f'''rearrangement <{barcodefile}.countfile -ref_file 3<{barcodefile}.reference -u -1 -v -3 -s0 -2 -qv -3 | correct_micro_homology.py {extup_left} {extdown_left} {extup_right} NGG | tee -a {barcodefile}.alg | awk -v OFS="\t" -v barcode={barcodegroup} -v cut1={extup_left} -v cut2={len(ref_left)+extup_right} 'NR%3==1{{print barcode, $0, cut1, cut2}}' >> {barcodefile}.table; echo >> {barcodefile}.alg''', shell=True, check=True)
    os.remove(f"{barcodefile}.countfile")
    os.remove(f"{barcodefile}.reference")
