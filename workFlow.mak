# Input: fastqFiles, spliterIndices, minScores, genome, bowtie2index, refFile, correctFile, ext1up, ext1down, ext2up, ext2down, s0, s1, s2, u, v, ru, rv, qu, qv, minToMapShear

comma = ,

spliterIndexBts = $(foreach dir,$(subst $(comma), ,$(spliterIndices)),$(addprefix $(dir)., 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2))

ifneq ($(targetSpliterFile),)
targetSpliterIndex = $(addprefix $(targetSpliterFile).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
endif
ifneq ($(pairSpliterFile),)
pairSpliterIndex = $(addprefix $(pairSpliterFile).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
endif

# outputDir includes the tail /
outputDir = $(dir $(word 1, $(subst $(comma), ,$(fastqFiles))))

$(outputDir)rearr.noDup: $(subst $(comma), ,$(fastqFiles))
	removeDuplicates.sh $^ >$@

%.1.bt2 %.2.bt2 %.3.bt2 %.4.bt2 %.rev.1.bt2 %.rev.2.bt2: %
	bowtie2-build $< $<

%.demultiplex: %.noDup $(spliterIndexBts)
	spliterIndices=$(spliterIndices) minScores=$(minScores) demultiplex.sh $< >$@

%.alg: %.post $(refFile) $(correctFile)
	rearrangement <$< 3<$(refFile) -s0 $(s0) -s1 $(s1) -s2 $(s2) -u $(u) -v $(v) -ru $(ru) -rv $(rv) -qu $(qu) -qv $(qv) | gawk -f correct_micro_homology.awk -- $(refFile) $(correctFile) >$@

# the followings are specific to sx data
%.ref: % $(genome) $(addprefix $(bowtie2index).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
	getSxCsvFileRef.sh $< $(genome) $(bowtie2index) $(ext1up) $(ext1down) $(ext2up) $(ext2down) >$@

%.correct: %.ref
	yes up | head -n$(shell wc -l <$<) >$@

%.post: %.demultiplex
	sxCutR2AdapterFilterCumulate.sh $< $(minToMapShear) >$@

%.target.fa %.pair.fa: %
	sxExtractSpliter.sh $< >$<.target.fa 3>$<.pair.fa

.PRECIOUS: $(outputDir)rearr.noDup %.1.bt2 %.demultiplex %.alg %.ref %.correct %.post %.csv.target.fa
