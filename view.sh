#!/bin/bash
# usage: view.sh correctfile readnum
# display the read alignments for the first readnum (default: 50) reads in ANSI fomrat

correctfile=$1 # alignment correct by micro homology (input.correct.cut.ext1.ext2)
read cut ext1 ext2 <<<$(echo $correctfile | awk -F "." -v OFS=" " '{print $(NF-2), $(NF-1), $NF}')
readnum=${2:-50} # read number to display (default: 50)
head -n$(($readnum * 3)) <$correctfile | align_align.py $(($cut + $ext1)) $cut $(($cut + $ext1 + $ext2)) | less -SNR # view the first 50 (150/3) read alignments