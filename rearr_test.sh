#!/bin/bash
project_path="$(dirname $(realpath $0))"
ref1=$(sed -n '2p' $project_path/test/ref12.fa)
ref2=$(sed -n '4p' $project_path/test/ref12.fa)
$project_path/rearr_run.sh $project_path/test/random.fq $ref1 $ref2 100 100 NGG NGG
$project_path/rearr_render.sh $project_path/test/random.fq.alg.100.100.200
$project_path/rearr_view.sh $project_path/test/random.fq.alg.100.100.200