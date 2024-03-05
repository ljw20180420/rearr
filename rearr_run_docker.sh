#!/bin/bash
#usage: rearr_run_docker.sh input ref1 ref2 cut1 cut2 NGGCCNtype1 NGGCCNtype2
# input can be (compressed) fasta/fastq

input=$1 # the input file
ref1=$2 # reference1
ref2=$3 # reference2
cut1=$4 # cut point of reference1
cut2=$5 # cut point of reference2
NGGCCNtype1=$6 # NGGCCNtype of reference1
NGGCCNtype2=$7 # NGGCCNtype of reference2
if [ ! -s $input ]
then
    echo "$input does not exist or is empty"
    exit 1
fi
project_path="$(dirname $(realpath $0))"
inputpath="$(dirname $(realpath $input))"
docker run --rm --mount type=bind,src=${inputpath},dst="/app/data" rearr:auto ./rearr_run.sh ./data/$(basename $input) $ref1 $ref2 $cut1 $cut2 $NGGCCNtype1 $NGGCCNtype2