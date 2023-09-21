#!/home/xiaoge/miniconda3/envs/crispr/bin/python
import sys, os, subprocess, Bio.Seq
barcodefile = sys.argv[1]
scaffold = Bio.Seq.Seq("gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc").reverse_complement().__str__().upper()
if not ".fq." in barcodefile and not ".fastq." in barcodefile:
    raise Exception("barcode file must contain fq or fastq")
csvpos = barcodefile.find(".fq.") + 4 if ".fq." in barcodefile else barcodefile.find(".fastq.") + 7
csvfile = os.path.join(os.path.dirname(barcodefile), barcodefile[csvpos:-8])
extup_left, extdown_left, extup_right, extdown_right = 50, 50, 50, 100
with os.popen(f'''cut -f1-2,6,8 {barcodefile}''', "r") as bcd, os.popen(f'''perl -anF, -E 'if($.>1){{$F[11]=~tr/ACGT/TGCA/; say $F[9] , "\t" , scalar reverse $F[11]}}' {csvfile} | sort -k2,2 | cut -f1 | bowtie2 -x /home/ljw/hg19_with_bowtie2_index/hg19 -r -U - 2> /dev/null | samtools view''', "r") as bf, os.popen(f'''tail -n+2 {csvfile} | cut -d',' -f12 | tr ACGT TGCA | rev | sort ''', "r") as rcb:
    _ = subprocess.run(f">{barcodefile}.alg", shell=True, check=True)
    bcdline = bcd.readline()
    if not bcdline:
        raise Exception("*.fq/fastq.barcode is empty")
    count, seq, barcodegroup, barcode_end = bcdline.rstrip().split()
    barcode_end = int(barcode_end)
    for line in bf:
        barcode = rcb.readline().rstrip()
        if barcode != barcodegroup:
            continue
        cf = open(f"barcode/countfile", "w")
        _ = cf.write(f"{seq[barcode_end + 3:]}\t{count}\n")
        _ = subprocess.run(f'''echo "{barcodegroup}" >> {barcodefile}.alg''', shell=True, check=True)
        qname, flag, chr, pos, _, CIGAR, _ = line.split("\t", 6)
        flag, pos, tglen = int(flag), int(pos), int(CIGAR[:-1])
        tgstrand = "-" if (flag // 16) % 2 else "+"
        cut = pos - 1 + 16 + 6 if tgstrand == "+" else pos - 1 + tglen - 16 - 6
        rfstrand = "+" if tgstrand == "-" else "-" # the strand of target is constract to that of barcode and reference
        if rfstrand == "-":
            eul, edl, eur, edr = extdown_left, extup_left, extdown_right, extup_right
        else:
            eul, edl, eur, edr = extup_left, extdown_left, extup_right, extdown_right
        ref_left = subprocess.run(f'''bedtools getfasta -bedOut -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed <(echo "{chr}"$'\t'"{cut-eul}"$'\t'"{cut+edl}"$'\t'{qname}$'\t'"."$'\t'"{rfstrand}")''', shell=True, check=True, capture_output=True, executable="/bin/bash").stdout.decode().split()[-1].rstrip()
        ref_right = subprocess.run(f'''bedtools getfasta -bedOut -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed <(echo "{chr}"$'\t'"{cut-eur}"$'\t'"{cut+edr}"$'\t'{qname}$'\t'"."$'\t'"{rfstrand}")''', shell=True, check=True, capture_output=True, executable="/bin/bash").stdout.decode().split()[-1].rstrip()
        with open("barcode/reference", "w") as rd:
            _ = rd.write(f"0\t0\t0\n{ref_left}\n{extup_left}\t{extup_left}\t{extup_left}\n{extup_right}\t{extup_right}\t{extup_right}\n{ref_right}\n{extup_right+extdown_right}\t{extup_right+extdown_right}\t{extup_right+extdown_right}\n0\t0\t0\n{scaffold}\n{len(scaffold)}\t{len(scaffold)}\t{len(scaffold)}\n")
        
        while True:
            bcdline = bcd.readline()
            if not bcdline:
                break
            count, seq, barcode, barcode_end = bcdline.rstrip().split()
            barcode_end = int(barcode_end)
            if barcode != barcodegroup:
                barcodegroup = barcode
                break
            _ = cf.write(f"{seq[barcode_end + 3:]}\t{count}\n")
        cf.close()
        _ = subprocess.run(f'''Rearrangement/build/rearrangement -file barcode/countfile -ref_file barcode/reference -ALIGN_MAX 1 -THR_NUM 1 >> {barcodefile}.alg; echo >> {barcodefile}.alg''', shell=True, check=True) # if thread number (-THR_NUM) set larger than 1, then the output does not keep order
    cf.close()
    _ = subprocess.run(f'''Rearrangement/build/rearrangement -file barcode/countfile -ref_file barcode/reference -ALIGN_MAX 1 -THR_NUM 1 >> {barcodefile}.alg; echo >> {barcodefile}.alg''', shell=True, check=True)