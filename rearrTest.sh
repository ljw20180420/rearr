#!/bin/bash

# Usage: ./rearrTest.sh
# To replace args (say makeTarget): makeTarget=test/A2-g1n-3.R2.fq.count ./rearrTest.sh [test]

if [[ $1 == "test" ]]
then
    targetSpliterFile=${targetSpliterFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.target.fa}
    pairSpliterFile=${pairSpliterFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.pair.fa}
fi

# The following parameters should be replaced.
makeTarget=${makeTarget:-test/A2-g1n-3.R2.fq.gz.alg}
pairFile=${pairFile:-test/A2-g1n-3.fq.gz}
minScoreTarget=${minScoreTarget:-30}
minScorePair=${minScorePair:-100}
minToMapShear=${minToMapShear:-30}
refFile=${refFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref}
ext1up=${ext1up:-50}
ext1down=${ext1down:-0}
ext2up=${ext2up:-10}
ext2down=${ext2down:-100}

# The following parameters are default in most cases.
genome=${genome:-genome/genome.fa}
bowtie2index=${bowtie2index:-genome/genome}
s0=${s0:--6}
s1=${s1:-4}
s2=${s2:-2}
u=${u:--3}
v=${v:--9}
ru=${ru:-0}
rv=${rv:-0}
qu=${qu:-0}
qv=${qv:--5}

make -f workFlow.mak $makeTarget pairFile=$pairFile targetSpliterFile=$targetSpliterFile pairSpliterFile=$pairSpliterFile genome=$genome bowtie2index=$bowtie2index refFile=$refFile minScoreTarget=$minScoreTarget minScorePair=$minScorePair s0=$s0 s1=$s1 s2=$s2 u=$u v=$v ru=$ru rv=$rv qu=$qu qv=$qv minToMapShear=$minToMapShear
