#!/bin/bash
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" mode reflen probability readnum software usertime systime realtime memory query refup refdown queryup querydown >bench/benchresult
for mode in single double
do
    for reflen in 100 200 300 400 500
    do
        for probability in 0.01 0.02 0.03 0.04 0.05
        do
            bench/tools/bench.sh $mode $reflen $probability 3>>bench/benchresult
        done
    done
done