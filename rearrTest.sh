#!/bin/bash

# Usage: ./test.sh
# To replace args (say makeTarget): makeTarget=test/A2-g1n-3.R2.fq.count ./test.sh

# The following parameters should be replaced.
makeTarget=${makeTarget:-test/A2-g1n-3.R2.fq.alg}
pairFile=${pairFile:-test/A2-g1n-3.fq}
targetSpliterFile=${targetSpliterFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.primer+barcode.fa}
pairSpliterFile=${pairSpliterFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.adapter+sgRNA+scaffold.fa}
refFile=${refFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref}

# The following parameters are default in most cases.
genome=${genome:-hg19/hg19.fa}
bowtie2index=${bowtie2index:-hg19/hg19}
minScoreTarget=${minScoreTarget:-30}
minScorePair=${minScorePair:-100}
s0=${s0:--6}
s1=${s1:-4}
s2=${s2:-2}
u=${u:--3}
v=${v:--9}
ru=${ru:-0}
rv=${rv:-0}
qu=${qu:-0}
qv=${qv:--5}
minToMapShear=${minToMapShear:-30}

make -f workFlow.mak $makeTarget pairFile=$pairFile targetSpliterFile=$targetSpliterFile pairSpliterFile=$pairSpliterFile genome=$genome bowtie2index=$bowtie2index refFile=$refFile minScoreTarget=$minScoreTarget minScorePair=$minScorePair s0=$s0 s1=$s1 s2=$s2 u=$u v=$v ru=$ru rv=$rv qu=$qu qv=$qv minToMapShear=$minToMapShear
