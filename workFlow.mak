# Input: fastqFiles, markerIndices, minScores, refFile, directionFile, ext1up, ext1down, ext2up, ext2down, s0, s1, s2, u, v, ru, rv, qu, qv

comma = ,

markerIndexBts = $(foreach dir,$(subst $(comma), ,$(markerIndices)),$(addprefix $(dir)., 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2))

# outputDir includes the tail /
outputDir = $(dir $(word 1, $(subst $(comma), ,$(fastqFiles))))

RN = $(length $(subst $(comma), ,$(fastqFiles)))

$(outputDir)rearr.noDup: $(subst $(comma), ,$(fastqFiles))
	removeDuplicates.sh $^ >$@

%.1.bt2 %.2.bt2 %.3.bt2 %.4.bt2 %.rev.1.bt2 %.rev.2.bt2: %
	bowtie2-build -q $< $<

%.demultiplex: %.noDup $(markerIndexBts)
	markerIndices=$(markerIndices) minScores=$(minScores) demultiplex.sh $< >$@

%.alg: %.post $(refFile) $(directionFile)
	rearrangement <$< 3<$(refFile) -s0 $(s0) -s1 $(s1) -s2 $(s2) -u $(u) -v $(v) -ru $(ru) -rv $(rv) -qu $(qu) -qv $(qv) | gawk -f correct_micro_homology.awk -- $(refFile) $(directionFile) >$@

%.direct: %
	gawk '{for (i = 1; i < NF / 3 - 1; ++i) printf("up\t"); printf("up\n");}' $< >$@

%.post: %.demultiplex
	cut -f1,$(shell expr $(RN) + 1)-$(shell expr $(RN) + 2) $< >$@

.PRECIOUS: $(outputDir)rearr.noDup %.1.bt2 %.demultiplex %.alg %.ref %.direct %.post