# Input: pairFile, targetSpliterFile, pairSpliterFile, genome, bowtie2index, refFile, minScoreTarget, minScorePair, s0, s1, s2, u, v, ru, rv, qu, qv
# sxInput: minToMapShear

%.count: % $(pairFile)
	removeDuplicates.sh $< $(pairFile) >$@

%.1.bt2 %.2.bt2 %.3.bt2 %.4.bt2 %.rev.1.bt2 %.rev.2.bt2: %
	bowtie2-build $< $<

%.demultiplex: %.count $(addprefix $(targetSpliterFile).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2) $(addprefix $(pairSpliterFile).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
	demultiplex.sh $< $(targetSpliterFile) $(pairSpliterFile) $(minScoreTarget) $(minScorePair) >$@

%.alg: %.post $(refFile)
	rearrangement <$< 3<$(refFile) -s0 $(s0) -s1 $(s1) -s2 $(s2) -u $(u) -v $(v) -ru $(ru) -rv $(rv) -qu $(qu) -qv $(qv) | gawk -f correct_micro_homology.awk -- $(refFile) NGG NGG >$@

# the followings are specific to sx data
%.ref: % $(genome) $(addprefix $(bowtie2index).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
	getSxCsvFileRef.sh $< $(genome) $(bowtie2index) >$@

%.post: %.demultiplex
	sxCutR2AdapterFilterCumulate.sh $< $(minToMapShear) >$@

%.csv.primer+barcode.fa %.csv.adapter+sgRNA+scaffold.fa: %.csv
	sxExtractSpliter.sh $<

.PRECIOUS: %.count %.1.bt2 %.demultiplex %.alg %.ref %.post %.csv.primer+barcode.fa
