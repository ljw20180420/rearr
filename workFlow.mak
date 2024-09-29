# Input: pairFile, targetSpliterFile, pairSpliterFile, genome, bowtie2index, refFile, ext1up, ext1down, ext2up, ext2down, minScoreTarget, minScorePair, s0, s1, s2, u, v, ru, rv, qu, qv, minToMapShear

ifneq ($(targetSpliterFile),)
targetSpliterIndex = $(addprefix $(targetSpliterFile).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
endif
ifneq ($(pairSpliterFile),)
pairSpliterIndex = $(addprefix $(pairSpliterFile).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
endif

%.count: % $(pairFile)
	removeDuplicates.sh $< $(pairFile) >$@

%.1.bt2 %.2.bt2 %.3.bt2 %.4.bt2 %.rev.1.bt2 %.rev.2.bt2: %
	bowtie2-build $< $<

%.demultiplex: %.count $(targetSpliterIndex) $(pairSpliterIndex)
	spliterTarget=$(targetSpliterFile) spliterPair=$(pairSpliterFile) minScoreTarget=$(minScoreTarget) minScorePair=$(minScorePair) demultiplex.sh $< >$@

%.alg: %.post $(refFile)
	rearrangement <$< 3<$(refFile) -s0 $(s0) -s1 $(s1) -s2 $(s2) -u $(u) -v $(v) -ru $(ru) -rv $(rv) -qu $(qu) -qv $(qv) | gawk -f correct_micro_homology.awk -- $(refFile) NGG NGG >$@

# the followings are specific to sx data
%.ref: % $(genome) $(addprefix $(bowtie2index).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
	getSxCsvFileRef.sh $< $(genome) $(bowtie2index) $(ext1up) $(ext1down) $(ext2up) $(ext2down) >$@

%.post: %.demultiplex
	sxCutR2AdapterFilterCumulate.sh $< $(minToMapShear) >$@

%.target.fa %.pair.fa: %
	sxExtractSpliter.sh $< >$<.target.fa 3>$<.pair.fa

.PRECIOUS: %.count %.1.bt2 %.demultiplex %.alg %.ref %.post %.csv.target.fa
