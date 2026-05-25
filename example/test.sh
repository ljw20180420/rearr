#!/bin/bash

# 切换运行路径到脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

makeTarget="example.alg" \
fastqFiles="example.fq" \
markerIndices="example.fa" \
minScores=20 \
refFile="example.ref" \
runWorkFlow.sh
