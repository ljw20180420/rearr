#!/bin/bash

# Change directory to the script directory.
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p complexity

ref_num=10
query_per_ref=1000

printf "reference length\tprobability\tpreprocess time\ttotal time\talign time\tDP time\tbacktrack time\treal time\tuser time\tsystem time\tmaximum memory\n" >complexity.tsv

for ref_len in 32 64 128 256 512 1024 2048 4096
do
    for probability in 0.01 0.02 0.03 0.04 0.05
    do
        utils/random_seq_methods.py $ref_len $ref_num $query_per_ref $probability double >complexity/ref.ref 3>complexity/query.post 4>/dev/null

        # time is a shell keyword, explicitly specify /bin/time to use the time executable
        /bin/time -f '%e\t%U\t%S\t%M' -o complexity/time.linux ../core/Rearrangement/build/rearrangement <complexity/query.post 3<complexity/ref.ref >/dev/null 2>complexity/time.tab

        printf "%d\t%f\t" $ref_len $probability >>complexity.tsv
        <complexity/time.tab sed -r 's/^[^\.0-9]+([\.0-9]+)s$/\1/' | sed -z 's/\n/\t/g' >>complexity.tsv
        cat complexity/time.linux >>complexity.tsv
    done
done

rm -rf complexity
