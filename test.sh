#!/bin/bash
targetSpliterFile=sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.primer+barcode.fa
pairSpliterFile=sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.adapter+sgRNA+scaffold.fa
refFile=sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref
make -f workFlow.mak test/A2-g1n-3.R2.fq.alg pairFile=test/A2-g1n-3.fq targetSpliterFile=${targetSpliterFile} pairSpliterFile=${pairSpliterFile} genomeRef=hg19/hg19.fa refFile=${refFile} minScoreTarget=30 minScorePair=100 s0=-6 s1=4 s2=2 u=-3 v=-9 ru=0 rv=0 qu=0 qv=-5 minToMapShear=30
