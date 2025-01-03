# Input: fastqFiles, spliterIndices, minScores, genome, bowtie2index, refFile, correctFile, ext1up, ext1down, ext2up, ext2down, s0, s1, s2, u, v, ru, rv, qu, qv, minToMapShear

comma = ,

spliterIndexBts = $(foreach dir,$(subst $(comma), ,$(spliterIndices)),$(addprefix $(dir)., 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2))

# outputDir includes the tail /
outputDir = $(dir $(word 1, $(subst $(comma), ,$(fastqFiles))))

$(outputDir)rearr.noDup: $(subst $(comma), ,$(fastqFiles))
	removeDuplicates.md $^ >$@

%.1.bt2 %.2.bt2 %.3.bt2 %.4.bt2 %.rev.1.bt2 %.rev.2.bt2: %
	bowtie2-build $< $<

%.demultiplex: %.noDup $(spliterIndexBts)
	spliterIndices=$(spliterIndices) minScores=$(minScores) demultiplex.md $< >$@

%.alg: %.post $(refFile) $(correctFile)
	rearrangement <$< 3<$(refFile) -s0 $(s0) -s1 $(s1) -s2 $(s2) -u $(u) -v $(v) -ru $(ru) -rv $(rv) -qu $(qu) -qv $(qv) | gawk -f correct_micro_homology.awk -- $(refFile) $(correctFile) >$@

# the followings are specific to sx data
%.ref: % $(genome) $(addprefix $(bowtie2index).,1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)
	getSxCsvFileRef.md $< $(genome) $(bowtie2index) $(ext1up) $(ext1down) $(ext2up) $(ext2down) >$@

%.correct: %.ref
	gawk '{for (i = 1; i < NF / 3 - 1; ++i) printf("up\t"); printf("up\n");}' $< >$@

%.post: %.demultiplex
	sxCutR2AdapterFilterCumulate.md $< $(minToMapShear) >$@

%.target.fa %.pair.fa: %
	sxExtractSpliter.md $< >$<.target.fa 3>$<.pair.fa

.PRECIOUS: $(outputDir)rearr.noDup %.1.bt2 %.demultiplex %.alg %.ref %.correct %.post %.target.fa %.pair.fa
