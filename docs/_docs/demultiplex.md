---
title: "Demultiplex"
permalink: /docs/demultiplex/
toc: true
---

# Usage: spliterTarget=spliterTarget spliterPair=spliterPair minScoreTarget=minScoreTarget minScorePair=minScorePair demultiplex.sh inputFile >demultiplexFile
# spliter(Target|Pair) is a bowtie2 index used to split (Target|pair)
# minScore(Target|Pair) thres the match between (Target|pair) and spliter(Target|Pair)

# inputFile format
# Target|pair|count
# Target is the read to map (either R1 or R2)
# pair is the read paired with Target
# count is the duplicate number of (Target, pair)

The `Nth` sequence in an `stdout` line, say `seqN`, is