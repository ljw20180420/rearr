#!/bin/bash

javac -g \
	-d paper/tools/SIQ-1.8/SIQ/bin \
	paper/tools/SIQ-1.8/SIQ/src/my/my_SIQ.java
jdb -classpath .:/home/ljw/new_fold/projects/rearr/paper/tools/SIQ-1.8/SIQ/libraries/* my.my_SIQ /tmp/tmp2wpd8u68/query.fq /tmp/tmp2wpd8u68/ref.fa TCGGAAGGCTAGGAG TTTGCATTCCTGCCG /tmp/tmp2wpd8u68/output.tsv
