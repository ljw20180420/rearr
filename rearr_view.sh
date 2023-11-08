#!/bin/bash
# usage: view.sh algfile [readnum]
# example: view.sh ljhlyz/AN1-SG4-M1B-1-1_R1.fq.gz.alg.75.30.30 50

# display the first readnum (default: 50) read's alignments in ANSI fomrat, the cut points and random insertions are aligned

algfile=$1 # alignments (input.alg.cut.ext1.ext2)
read cut ext1 ext2 <<<$(echo $algfile | awk -F "." -v OFS=" " '{print $(NF-2), $(NF-1), $NF}')
readnum=${2:-50} # read number to display (default: 50)
head -n$(($readnum * 3)) <$algfile | align_align.py $(($cut + $ext1)) $cut $(($cut + $ext1 + $ext2)) | less-643/less -SNR --header 1 --no-number-headers # view the first 50 (150/3) read alignments