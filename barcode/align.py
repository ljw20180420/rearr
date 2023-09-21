import Bio.Seq, os, subprocess

barcodefile = "barcode/A2_TEST.fq.barcode"
ext = 1000
csvfile = "barcode/final_hgsgrna_libb_all_0811-NGG.csv"
os.makedirs("output", exist_ok=True)


with os.popen(f'''cut -f1-2,6 {barcodefile}''', "r") as bcd, os.popen(f'''perl -anF, -E 'if($.>1){{$F[11]=~tr/ACGT/TGCA/; say $F[9] , "\t" , scalar reverse $F[11]}}' {csvfile} | sort -k2,2 | cut -f1 | bowtie2 -x /home/ljw/hg19_with_bowtie2_index/hg19 -r -U - 2> /dev/null | samtools view''', "r") as bf, os.popen(f'''tail -n+2 {csvfile} | cut -d',' -f12 | tr ACGT TGCA | rev | sort ''', "r") as rcb:
    bcdline = bcd.readline()
    count, seq, barcodegroup = bcdline.rstrip().split() 
    for line in bf:
        if not bcdline:
            break
        cf = os.popen(f"sort -k1,1nr > output/countfile.{barcodegroup}", "w")
        _ = cf.write(f"{seq}\t{count}\n")
        barcode = rcb.readline().rstrip()
        if barcode != barcodegroup:
            continue
        qname, flag, chr, pos, _, CIGAR, _ = line.split("\t", 6)
        flag, pos, tglen = int(flag), int(pos), int(CIGAR[:-1])
        tgstrand = "-" if (flag // 16) % 2 else "+"
        cut = pos + 16 + 6 if tgstrand == "+" else pos + tglen - 16 - 6
        rfstrand = "+" if tgstrand == "-" else "-" # the strand of target is constract to that of barcode and reference
        ref = next(subprocess.Popen(f'''bedtools getfasta -bedOut -s -fi /home/ljw/hg19_with_bowtie2_index/hg19.fa -bed <(echo "{chr}"$'\t'"{cut-ext}"$'\t'"{cut+ext}"$'\t'{qname}$'\t'"."$'\t'"{rfstrand}")''', shell=True, stdout=subprocess.PIPE, executable="/bin/bash").stdout).decode().split()[-1].rstrip()
        with open("output/reference", "w") as rd:
            _ = rd.write(f"0\t0\t0\n{ref}\n{ext}\t{ext}\t{ext}\n{ext}\t{ext}\t{ext}\n{ref}\n{2*ext}\t{2*ext}\t{2*ext}\n")
        
        while True:
            bcdline = bcd.readline()
            if not bcdline:
                break
            count, seq, barcode = bcdline.rstrip().split() 
            if barcode != barcodegroup:
                barcodegroup = barcode
                break
            _ = cf.write(f"{seq}\t{count}\n")
        cf.close()
        subprocess.check_output(f'''cd output; /home/ljw/wuqiang/lierlib/Rearrangement/build/rearrangement -file countfile.{barcodegroup} -ref_file reference -ALIGN_MAX 1 -THR_NUM 1''', shell=True)
    cf.close()
    subprocess.check_output(f'''cd output; /home/ljw/wuqiang/lierlib/Rearrangement/build/rearrangement -file countfile.{barcodegroup} -ref_file reference -ALIGN_MAX 1 -THR_NUM 1''', shell=True)