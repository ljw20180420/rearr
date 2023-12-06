#!/bin/bash
generate_random_DNA()
{
    # usage: generate_random_DNA length
    tr -dc ACGT </dev/urandom | head -c $1; echo ''
}


reflen=200
sglen=20
cut=$(expr $reflen / 2)
sgstart=$(expr $cut - $sglen + 3)
ext=$(expr $reflen / 20)
ref=$(generate_random_DNA $reflen)
ref=$(printf $ref | head -c$(($sgstart+$sglen+1)))GG$(printf $ref | tail -c$(($reflen-$sgstart-$sglen-3)))
sed '1i >ref' <<<$ref >test/ref.fa
sgRNA=$(printf $ref | cut -c$(($sgstart+1))-$(($sgstart+$sglen)))

readnum=10000
bench/tools/random_seq_methods.py $ref 0.02 $readnum >test/random.seq
cut -f1 test/random.seq | perl -nE 'say "\@seq" . $. . "\n" . $_ . "+\n" . "~" x (length($_) - 1)' > test/random.fq
cut -f2- test/random.seq > test/random.seq2; mv test/random.seq2 test/random.seq

rm -rf test/rearr; mkdir -p test/rearr
cp test/random.fq test/rearr/random.fq
rearr_run.sh test/rearr/random.fq $ref $sgRNA $ext $ext
tail -n+2 test/rearr/random.fq.table.$cut.$ext.$ext | awk -F "\t" -v OFS="\t" -v cut=$cut -v ext=$ext '{$10 -= 2 * ext; print $7,$10,$8, $11}' | paste <(cut -f1 test/rearr/random.fq.count) - | sort -k1,1 | join -t $'\t' -1 2 -2 1 <(sed -nr '1~4{s/^@//; N; s/\n/\t/; p}' test/rearr/random.fq | sort -k2,2) - | cut -f2- | sort -k1,1V > test/rearr.arr

join -t $'\t' -1 1 -2 1 <(sort -k1,1 test/random.seq) <(sort -k1,1 test/rearr.arr) | sort -k1,1V | bench/tools/compare_indel.py $ref >test/rearr.diff

bench/tools/draw_bench.py $readnum "test/rearr.diff"