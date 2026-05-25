#!/bin/bash

# 切换运行路径到脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

. lib.sh

# generate_seq_groupby_marker 3> test.fa 4> test.ref | make_fastq > test.fq

makeTarget="test.alg" \
fastqFiles="test.fq" \
markerIndices="test.fa" \
minScores=20 \
refFile="test.ref" \
runWorkFlow.sh
